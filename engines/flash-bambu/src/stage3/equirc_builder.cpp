#include <Rcpp.h>
#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <numeric>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "common/flat_csr_index.hpp"
#include "stage3/csr_splice_dist.hpp"

namespace bambu_next {
namespace {

struct SigVecHash {
    size_t operator()(const std::vector<int32_t> &v) const noexcept {
        uint64_t h = 0x9e3779b97f4a7c15ULL ^ v.size();
        for (int32_t x : v) {
            uint64_t k = static_cast<uint64_t>(static_cast<uint32_t>(x));
            k = (k ^ (k >> 30)) * 0xbf58476d1ce4e5b9ULL;
            h ^= k + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return static_cast<size_t>(h);
    }
};

// Cached immutable minimal equivalence class index for an extended FlatCSRAnnotation
struct ExtMinEqCache {
    const FlatCSRAnnotation *csr_ptr = nullptr;
    int32_t n_tx = 0;
    size_t eq_size = 0;

    std::vector<int32_t> tx_len;                  // size n_tx
    std::vector<std::vector<int32_t>> gene_txs;   // per gene: 0-based tx indices
    // Unique positive signatures (all positive 1-based txids, sorted)
    std::vector<std::vector<int32_t>> pos_sigs;
    std::unordered_map<std::vector<int32_t>, int32_t, SigVecHash> pos_sig_to_id;
    std::vector<std::vector<int32_t>> pos_sig_genes; // per pos_sig_id: sorted unique 0-based gene_ids
    std::vector<std::vector<int32_t>> gene_pos_sigs; // per gene: sorted unique pos_sig_id
};

static std::mutex g_ext_mineq_mtx;
static std::deque<std::shared_ptr<ExtMinEqCache>> g_ext_mineq_lru;

static std::shared_ptr<ExtMinEqCache> get_or_build_ext_mineq_cache(const FlatCSRAnnotation &ann) {
    {
        std::lock_guard<std::mutex> lk(g_ext_mineq_mtx);
        for (const auto &e : g_ext_mineq_lru) {
            if (e->csr_ptr == &ann && e->n_tx == ann.n_tx && e->eq_size == ann.eq_tx.size()) {
                return e;
            }
        }
    }

    auto cache = std::make_shared<ExtMinEqCache>();
    cache->csr_ptr = &ann;
    cache->n_tx = ann.n_tx;
    cache->eq_size = ann.eq_tx.size();
    cache->tx_len.assign(ann.n_tx, 0);
    cache->gene_txs.resize(ann.n_gene);
    cache->gene_pos_sigs.resize(ann.n_gene);
    cache->pos_sig_to_id.reserve(static_cast<size_t>(ann.n_tx) + 64);

    std::vector<int32_t> tmp_sig;
    for (int32_t t = 0; t < ann.n_tx; ++t) {
        int32_t tlen = 0;
        for (int32_t k = ann.tx_exon_ptr[t]; k < ann.tx_exon_ptr[t + 1]; ++k) {
            tlen += (ann.exon_end[k] - ann.exon_start[k] + 1);
        }
        cache->tx_len[t] = tlen;
        int32_t gid = ann.tx_geneid[t];
        cache->gene_txs[gid].push_back(t);

        int32_t eb = ann.eq_ptr[t], ee = ann.eq_ptr[t + 1];
        tmp_sig.resize(ee - eb);
        for (int32_t k = eb; k < ee; ++k) tmp_sig[k - eb] = ann.eq_tx[k] + 1;
        // eq_tx is already sorted, positive 1-based is still sorted
        auto it = cache->pos_sig_to_id.find(tmp_sig);
        int32_t ps_id = -1;
        if (it == cache->pos_sig_to_id.end()) {
            ps_id = static_cast<int32_t>(cache->pos_sigs.size());
            cache->pos_sigs.push_back(tmp_sig);
            cache->pos_sig_to_id.emplace(cache->pos_sigs.back(), ps_id);
            cache->pos_sig_genes.emplace_back();
        } else {
            ps_id = it->second;
        }
        cache->pos_sig_genes[ps_id].push_back(gid);
        cache->gene_pos_sigs[gid].push_back(ps_id);
    }

    for (auto &v : cache->pos_sig_genes) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    for (auto &v : cache->gene_pos_sigs) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    {
        std::lock_guard<std::mutex> lk(g_ext_mineq_mtx);
        for (const auto &e : g_ext_mineq_lru) {
            if (e->csr_ptr == &ann && e->n_tx == ann.n_tx && e->eq_size == ann.eq_tx.size()) {
                return e;
            }
        }
        if (g_ext_mineq_lru.size() >= 2) g_ext_mineq_lru.pop_front();
        g_ext_mineq_lru.push_back(cache);
    }
    return cache;
}

struct EqObsEntry {
    std::vector<int32_t> sig; // signed 1-based txid (* (-1)^equal), sorted ascending
    int32_t gene_id;
    double nobs;
    double rc_width;
    bool has_neg; // any equal == TRUE in signature
    bool has_pos; // any equal == FALSE in signature
};

struct EqRow {
    int32_t eq_local_id;
    int32_t tx_idx;     // 0-based into ann
    double nobs;
    double rc_width;
    bool equal;
    double aval;
    double full_aval;
    double unique_aval;
};

struct GeneSparseBlock {
    int32_t gene_id;
    int32_t n_tx;
    int32_t n_eq;
    int64_t tr_dim;
    double gene_k;
    std::vector<int32_t> active_txs; // sorted global 0-based tx_idx
    std::vector<double> eq_y;        // size n_eq (nobs / gene_k)
    std::vector<double> eq_nobs;     // size n_eq
    struct Edge {
        int32_t tx_l;
        int32_t eq_l;
        double aval;
        double full_aval;
        double uniq_aval;
    };
    std::vector<Edge> edges;
};

} // namespace
} // namespace bambu_next

extern "C" SEXP bambu_next_quantify_cpp(
    SEXP csr_ptr_sexp,
    SEXP q_seq_names_sexp,
    SEXP q_tx_ex_counts_sexp,
    SEXP q_ex_seq_1b_sexp,
    SEXP q_ex_start_sexp,
    SEXP q_ex_end_sexp,
    SEXP q_ex_strand_sexp,
    SEXP q_ex_rank_sexp,
    SEXP q_ex_end_rank_sexp,
    SEXP q_read_count_sexp,
    SEXP q_rc_names_sexp,
    SEXP isore_params_sexp,
    SEXP em_params_sexp,
    SEXP return_dist_table_sexp,
    SEXP n_threads_sexp
) {
    BEGIN_RCPP
    if (Rf_isNull(csr_ptr_sexp)) Rcpp::stop("FlatCSRAnnotation pointer is null in Stage 3 quantify");
    Rcpp::XPtr<bambu_next::FlatCSRAnnotation> ann(csr_ptr_sexp);
    auto cache = bambu_next::get_or_build_ext_mineq_cache(*ann);

    Rcpp::CharacterVector q_snames_r(q_seq_names_sexp);
    Rcpp::IntegerVector q_cnts_r(q_tx_ex_counts_sexp);
    Rcpp::IntegerVector q_es_r(q_ex_seq_1b_sexp);
    Rcpp::IntegerVector q_s_r(q_ex_start_sexp);
    Rcpp::IntegerVector q_e_r(q_ex_end_sexp);
    Rcpp::IntegerVector q_st_r(q_ex_strand_sexp);
    Rcpp::IntegerVector q_rk_r(q_ex_rank_sexp);
    Rcpp::IntegerVector q_erk_r(q_ex_end_rank_sexp);
    Rcpp::NumericVector q_rc_r(q_read_count_sexp);
    Rcpp::CharacterVector q_names_r(q_rc_names_sexp);

    Rcpp::List iso_p(isore_params_sexp);
    int32_t max_dist = Rcpp::as<int32_t>(iso_p["min_exonDistance"]);
    int32_t ps_dist = Rcpp::as<int32_t>(iso_p["min_primarySecondaryDist"]);
    int32_t ps_dist_se = Rcpp::as<int32_t>(iso_p["min_primarySecondaryDistStartEnd2"]);

    Rcpp::List em_p(em_params_sexp);
    bool d_mode = Rcpp::as<bool>(em_p["degradationBias"]);
    int maxiter = Rcpp::as<int>(em_p["maxiter"]);
    double conv = Rcpp::as<double>(em_p["conv"]);
    double minvalue = Rcpp::as<double>(em_p["minvalue"]);
    int sig_digit = Rcpp::as<int>(em_p["sig_digit"]);

    bool return_dist_table = Rcpp::as<bool>(return_dist_table_sexp);
    int n_threads = std::max(1, Rf_asInteger(n_threads_sexp));

    int32_t n_q = static_cast<int32_t>(q_cnts_r.size());
    int32_t n_tx = ann->n_tx;
    int32_t n_gene = ann->n_gene;

    std::vector<std::string> q_snames(q_snames_r.size());
    for (R_xlen_t i = 0; i < q_snames_r.size(); ++i) q_snames[i] = Rcpp::as<std::string>(q_snames_r[i]);

    std::vector<int32_t> q_ptr(n_q + 1);
    std::vector<int32_t> rc_first_ex_w(n_q, 0);
    std::vector<int32_t> rc_tot_w(n_q, 0);
    std::vector<double> rc_read_cnt(n_q, 0.0);
    int32_t off = 0;
    for (int32_t i = 0; i < n_q; ++i) {
        q_ptr[i] = off;
        off += q_cnts_r[i];
        rc_read_cnt[i] = q_rc_r[i];
    }
    q_ptr[n_q] = off;

    std::vector<int32_t> q_es(off), q_s(off), q_e(off), q_rk(off), q_erk(off);
    std::vector<int8_t> q_st(off);
    for (int32_t i = 0; i < off; ++i) {
        q_es[i] = q_es_r[i] - 1;
        q_s[i] = q_s_r[i];
        q_e[i] = q_e_r[i];
        q_st[i] = static_cast<int8_t>(q_st_r[i]);
        q_rk[i] = q_rk_r[i];
        q_erk[i] = q_erk_r[i];
    }
    for (int32_t q = 0; q < n_q; ++q) {
        int32_t tw = 0, fw = 0;
        for (int32_t k = q_ptr[q]; k < q_ptr[q + 1]; ++k) {
            int32_t w = q_e[k] - q_s[k] + 1;
            tw += w;
            if (q_rk[k] == 1) fw = w;
        }
        rc_first_ex_w[q] = fw;
        rc_tot_w[q] = tw;
    }

    // Step 1: C++ Stage 3 calculateDistToAnnotation across all Q read classes
    std::vector<std::vector<bambu_next::DistHit>> per_q;
    bambu_next::calculate_dist_to_annotation_core(
        *ann, q_snames, q_ptr, q_es, q_s, q_e, q_st, q_rk, q_erk,
        max_dist, ps_dist, ps_dist_se, false, n_threads, per_q
    );

    // Step 2: Assign geneId per readClass matching addGeneIdsToReadClassTable
    //   geneCount[gid] = sum(readCount) across ALL distTable rows of gid (including duplicates within same RC)
    std::vector<double> gene_dist_sum(n_gene, 0.0);
    for (int32_t q = 0; q < n_q; ++q) {
        double rc = rc_read_cnt[q];
        for (const auto &h : per_q[q]) {
            gene_dist_sum[ann->tx_geneid[h.subject_idx]] += rc;
        }
    }

    std::vector<int32_t> rc_gene_id(n_q, -1);
    std::vector<uint8_t> rc_any_comp(n_q, 0);
    std::vector<double> incomp_by_gene(n_gene, 0.0);

    for (int32_t q = 0; q < n_q; ++q) {
        const auto &hits = per_q[q];
        if (hits.empty()) continue;
        bool any_c = false;
        int32_t best_g = -1;
        double best_g_sum = -1.0;
        for (const auto &h : hits) {
            if (h.compatible) any_c = true;
            int32_t g = ann->tx_geneid[h.subject_idx];
            double gs = gene_dist_sum[g];
            if (gs > best_g_sum) { // strictly > selects row_number() == 1 on ties in distTable order!
                best_g_sum = gs;
                best_g = g;
            }
        }
        rc_gene_id[q] = best_g;
        rc_any_comp[q] = any_c ? 1 : 0;
        if (!any_c && best_g >= 0) {
            incomp_by_gene[best_g] += rc_read_cnt[q];
        }
    }

    // Step 3: Build observed equivalence classes (genEquiRCsBasedOnObservedReads + getUniCountPerEquiRC)
    struct ObsRCInfo {
        int32_t fw;
        int32_t tw;
        double rc;
        int32_t gene_id;
        bool has_neg;
        bool has_pos;
    };
    struct SigBucket {
        std::vector<int32_t> sig;
        int32_t gene_id;
        std::vector<ObsRCInfo> items;
    };
    std::vector<SigBucket> obs_buckets;
    // Key: sig -> unique sig id (preserving first-appearance order across q)
    std::vector<std::vector<int32_t>> sig_keys;
    std::vector<std::vector<int32_t>> sig_to_buckets;
    std::unordered_map<std::vector<int32_t>, int32_t, bambu_next::SigVecHash> obs_sig_to_id;
    obs_buckets.reserve(n_q / 2 + 64);
    sig_keys.reserve(n_q / 2 + 64);
    sig_to_buckets.reserve(n_q / 2 + 64);
    obs_sig_to_id.reserve(n_q + 64);

    std::vector<int32_t> tmp_sig;
    for (int32_t q = 0; q < n_q; ++q) {
        if (!rc_any_comp[q] || rc_gene_id[q] < 0) continue;
        int32_t gid = rc_gene_id[q];
        tmp_sig.clear();
        bool h_neg = false, h_pos = false;
        for (const auto &h : per_q[q]) {
            if (ann->tx_geneid[h.subject_idx] == gid) {
                int32_t val = h.equal ? -(h.subject_idx + 1) : (h.subject_idx + 1);
                tmp_sig.push_back(val);
                if (val < 0) h_neg = true;
                else h_pos = true;
            }
        }
        if (tmp_sig.empty()) continue;
        std::sort(tmp_sig.begin(), tmp_sig.end());

        auto it = obs_sig_to_id.find(tmp_sig);
        int32_t sid_idx = -1;
        int32_t b_idx = -1;
        if (it != obs_sig_to_id.end()) {
            sid_idx = it->second;
            for (int32_t idx : sig_to_buckets[sid_idx]) {
                if (obs_buckets[idx].gene_id == gid) {
                    b_idx = idx;
                    break;
                }
            }
        } else {
            sid_idx = static_cast<int32_t>(sig_keys.size());
            sig_keys.push_back(tmp_sig);
            sig_to_buckets.emplace_back();
            obs_sig_to_id.emplace(sig_keys.back(), sid_idx);
        }
        if (b_idx < 0) {
            b_idx = static_cast<int32_t>(obs_buckets.size());
            obs_buckets.push_back(SigBucket{tmp_sig, gid, {}});
            sig_to_buckets[sid_idx].push_back(b_idx);
        }
        obs_buckets[b_idx].items.push_back(ObsRCInfo{
            rc_first_ex_w[q], rc_tot_w[q], rc_read_cnt[q], gid, h_neg, h_pos
        });
    }

    // Aggregate getUniCountPerEquiRC per observed (sig, gene_id) in exact first-appearance order,
    // AND replicate addEmptyRC() appending any shared unobserved minEquiRC genes afterwards
    // so simplifyNames() gene_sid order (and assignGroups() tie-breaks) match native R 100%.
    std::vector<bambu_next::EqObsEntry> obs_eqs;
    std::vector<bambu_next::EqObsEntry> shared_extra_eqs;
    obs_eqs.reserve(obs_buckets.size());
    std::vector<uint8_t> gene_active(n_gene, 0);
    std::vector<int32_t> active_genes;
    active_genes.reserve(n_gene);

    for (size_t si = 0; si < sig_keys.size(); ++si) {
        const auto &sig = sig_keys[si];
        const auto &b_list = sig_to_buckets[si];
        bool grp_any_eq = (!sig.empty() && sig[0] < 0);
        // Aggregate over unique (fw, tw, rc, gid, anyEqual) per sig
        std::vector<ObsRCInfo> all_items;
        for (int32_t bi : b_list) {
            all_items.insert(all_items.end(), obs_buckets[bi].items.begin(), obs_buckets[bi].items.end());
        }
        std::sort(all_items.begin(), all_items.end(), [](const ObsRCInfo &a, const ObsRCInfo &b) {
            if (a.fw != b.fw) return a.fw < b.fw;
            if (a.tw != b.tw) return a.tw < b.tw;
            if (a.rc != b.rc) return a.rc < b.rc;
            return a.gene_id < b.gene_id;
        });
        double nobs = 0.0;
        double rcw = 0.0;
        for (size_t i = 0; i < all_items.size(); ++i) {
            if (i > 0 &&
                all_items[i].fw == all_items[i - 1].fw &&
                all_items[i].tw == all_items[i - 1].tw &&
                all_items[i].rc == all_items[i - 1].rc &&
                all_items[i].gene_id == all_items[i - 1].gene_id) {
                continue; // Replicate dplyr::distinct() in getUniCountPerEquiRC!
            }
            nobs += all_items[i].rc;
            double w = grp_any_eq ? static_cast<double>(all_items[i].tw) : static_cast<double>(all_items[i].fw);
            if (w > rcw) rcw = w;
        }
        bool h_neg = (!sig.empty() && sig[0] < 0);
        bool h_pos = (!sig.empty() && sig.back() > 0);

        // Observed genes first in exact b_list (q-appearance) order
        std::vector<int32_t> obs_gids;
        obs_gids.reserve(b_list.size());
        for (int32_t bi : b_list) {
            int32_t g = obs_buckets[bi].gene_id;
            if (std::find(obs_gids.begin(), obs_gids.end(), g) == obs_gids.end()) {
                obs_gids.push_back(g);
                obs_eqs.push_back(bambu_next::EqObsEntry{sig, g, nobs, rcw, h_neg, h_pos});
                if (nobs > 0.0 && !gene_active[g]) {
                    gene_active[g] = 1;
                    active_genes.push_back(g);
                }
            }
        }

        // Extra shared minEquiRC genes (appended after eqClassCount in full_join)
        std::vector<int32_t> extra_gids;
        if (!h_neg) {
            auto pit = cache->pos_sig_to_id.find(sig);
            if (pit != cache->pos_sig_to_id.end()) {
                for (int32_t g2 : cache->pos_sig_genes[pit->second]) extra_gids.push_back(g2);
            }
        } else {
            int32_t cand_t = (-sig[0]) - 1;
            if (cand_t >= 0 && cand_t < n_tx) {
                int32_t eb = ann->eq_ptr[cand_t], ee = ann->eq_ptr[cand_t + 1];
                if (static_cast<size_t>(ee - eb) == sig.size()) {
                    std::vector<int32_t> full_s(ee - eb);
                    for (int32_t k = eb; k < ee; ++k) {
                        int32_t u = ann->eq_tx[k];
                        full_s[k - eb] = (u == cand_t) ? -(u + 1) : (u + 1);
                    }
                    std::sort(full_s.begin(), full_s.end());
                    if (full_s == sig) extra_gids.push_back(ann->tx_geneid[cand_t]);
                }
            }
        }
        std::sort(extra_gids.begin(), extra_gids.end());
        extra_gids.erase(std::unique(extra_gids.begin(), extra_gids.end()), extra_gids.end());
        for (int32_t g : extra_gids) {
            if (std::find(obs_gids.begin(), obs_gids.end(), g) == obs_gids.end()) {
                shared_extra_eqs.push_back(bambu_next::EqObsEntry{sig, g, nobs, rcw, h_neg, h_pos});
            }
        }
    }
    for (const auto &se : shared_extra_eqs) {
        obs_eqs.push_back(se);
        if (se.nobs > 0.0 && !gene_active[se.gene_id]) {
            gene_active[se.gene_id] = 1;
            active_genes.push_back(se.gene_id);
        }
    }

    // Group observed EqObsEntry by gene_id
    std::vector<std::vector<int32_t>> gene_obs_idx(n_gene);
    for (size_t i = 0; i < obs_eqs.size(); ++i) {
        gene_obs_idx[obs_eqs[i].gene_id].push_back(static_cast<int32_t>(i));
    }

    std::vector<double> g_nobs_tot(active_genes.size(), 0.0);
    std::vector<double> g_dobs_tot(active_genes.size(), 0.0);
    std::vector<int32_t> g_maxlen(active_genes.size(), 1);
    std::vector<std::vector<bambu_next::EqRow>> gene_rows(active_genes.size());
    std::vector<std::vector<double>> gene_eq_nobs(active_genes.size());
    std::vector<std::vector<double>> gene_eq_rcw(active_genes.size());

    // Build per-active-gene equivalence table in parallel
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(n_threads) if (active_genes.size() > 256 && n_threads > 1)
#endif
    for (size_t gi = 0; gi < active_genes.size(); ++gi) {
        int32_t gid = active_genes[gi];
        std::unordered_map<std::vector<int32_t>, int32_t, bambu_next::SigVecHash> local_sig_map;
        const auto &oidx = gene_obs_idx[gid];
        local_sig_map.reserve((oidx.size() + cache->gene_txs[gid].size() * 2 + 4) * 2);

        double tot_nobs = 0.0;
        double par_nobs = 0.0;
        int32_t max_len = 0;
        auto &rows = gene_rows[gi];
        auto &eq_nobs = gene_eq_nobs[gi];
        auto &eq_rcw = gene_eq_rcw[gi];
        int32_t next_eq = 0;

        auto ensure_sig = [&](const std::vector<int32_t> &sig, double nobs, double rcw) -> int32_t {
            auto it = local_sig_map.find(sig);
            if (it == local_sig_map.end()) {
                int32_t id = next_eq++;
                local_sig_map.emplace(sig, id);
                eq_nobs.push_back(nobs);
                eq_rcw.push_back(rcw);
                return id;
            }
            return it->second;
        };

        // Add observed eq classes for gid
        for (int32_t oi : oidx) {
            const auto &oe = obs_eqs[oi];
            int32_t eq_l = ensure_sig(oe.sig, oe.nobs, oe.rc_width);
            tot_nobs += oe.nobs;
            if (oe.has_pos) par_nobs += oe.nobs;
            for (int32_t stx : oe.sig) {
                bool eq_flag = (stx < 0);
                int32_t tx = std::abs(stx) - 1;
                if (cache->tx_len[tx] > max_len) max_len = cache->tx_len[tx];
                rows.push_back(bambu_next::EqRow{eq_l, tx, oe.nobs, oe.rc_width, eq_flag, 1.0, 0.0, 0.0});
            }
        }

        // Add unobserved positive minEquiRC signatures for gid
        std::vector<int32_t> sig_buf;
        for (int32_t ps_id : cache->gene_pos_sigs[gid]) {
            const auto &psig = cache->pos_sigs[ps_id];
            if (local_sig_map.find(psig) != local_sig_map.end()) continue; // already observed!
            int32_t eq_l = ensure_sig(psig, 0.0, 0.0);
            for (int32_t stx : psig) {
                int32_t tx = stx - 1;
                if (cache->tx_len[tx] > max_len) max_len = cache->tx_len[tx];
                rows.push_back(bambu_next::EqRow{eq_l, tx, 0.0, 0.0, false, 1.0, 0.0, 0.0});
            }
        }

        // Add unobserved full-length minEquiRCTemp signatures for each t in gene_txs[gid]
        for (int32_t t : cache->gene_txs[gid]) {
            int32_t eb = ann->eq_ptr[t], ee = ann->eq_ptr[t + 1];
            sig_buf.resize(ee - eb);
            for (int32_t k = eb; k < ee; ++k) {
                int32_t u = ann->eq_tx[k];
                sig_buf[k - eb] = (u == t) ? -(u + 1) : (u + 1);
            }
            std::sort(sig_buf.begin(), sig_buf.end());
            if (local_sig_map.find(sig_buf) != local_sig_map.end()) continue;
            int32_t eq_l = ensure_sig(sig_buf, 0.0, 0.0);
            for (int32_t stx : sig_buf) {
                bool eq_flag = (stx < 0);
                int32_t tx = std::abs(stx) - 1;
                if (cache->tx_len[tx] > max_len) max_len = cache->tx_len[tx];
                rows.push_back(bambu_next::EqRow{eq_l, tx, 0.0, 0.0, eq_flag, 1.0, 0.0, 0.0});
            }
        }

        g_nobs_tot[gi] = tot_nobs;
        g_dobs_tot[gi] = par_nobs;
        g_maxlen[gi] = std::max(1, max_len);
    }

    // Calculate degradation rate (d_rate) identically to calculateDegradationRate()
    double d_rate = NA_REAL;
    if (d_mode) {
        std::vector<double> filtered_rates;
        std::vector<double> all_rates;
        filtered_rates.reserve(active_genes.size());
        all_rates.reserve(active_genes.size());
        for (size_t gi = 0; gi < active_genes.size(); ++gi) {
            double nobs = g_nobs_tot[gi];
            double dobs = g_dobs_tot[gi];
            if (nobs <= 0.0) continue;
            double val = (dobs / nobs) * 1000.0 / static_cast<double>(g_maxlen[gi]);
            all_rates.push_back(val);
            if (nobs >= 30.0 && (nobs - dobs) >= 5.0) {
                filtered_rates.push_back(val);
            }
        }
        auto &use_rates = filtered_rates.empty() ? all_rates : filtered_rates;
        if (filtered_rates.empty()) {
            Rcpp::message(Rf_mkString("There is not enough read count and full length coverage!\n            Hence degradation rate is estimated using all data!"));
        }
        if (!use_rates.empty()) {
            std::sort(use_rates.begin(), use_rates.end());
            size_t m = use_rates.size();
            d_rate = (m % 2 == 1) ? use_rates[m / 2] : 0.5 * (use_rates[m / 2 - 1] + use_rates[m / 2]);
        }
    }

    // Step 5: Compute aval, filter to sum(nobs) > 0 per (gene_sid, txid), build sparse blocks
    std::vector<bambu_next::GeneSparseBlock> g_blocks(active_genes.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(n_threads) if (active_genes.size() > 256 && n_threads > 1)
#endif
    for (size_t gi = 0; gi < active_genes.size(); ++gi) {
        auto &rows = gene_rows[gi];
        int32_t n_eq = static_cast<int32_t>(gene_eq_nobs[gi].size());

        // Count unique tx per eq_l for multi_align, and track sum(nobs) per tx & per eq
        std::sort(rows.begin(), rows.end(), [](const bambu_next::EqRow &a, const bambu_next::EqRow &b) {
            if (a.tx_idx != b.tx_idx) return a.tx_idx < b.tx_idx;
            return a.eq_local_id < b.eq_local_id;
        });

        std::vector<int32_t> eq_deg(n_eq, 0);
        std::vector<uint8_t> eq_all_not_eq(n_eq, 1);
        for (const auto &r : rows) {
            eq_deg[r.eq_local_id]++;
            if (r.equal) eq_all_not_eq[r.eq_local_id] = 0;
        }

        // Group by tx_idx -> compute sum_not_eq_rcw and sum_nobs
        std::vector<int32_t> active_txs;
        std::unordered_map<int32_t, int32_t> tx_to_local;
        size_t i = 0;
        while (i < rows.size()) {
            size_t j = i;
            int32_t tx = rows[i].tx_idx;
            double tx_nobs = 0.0;
            double sum_noteq_rcw = 0.0;
            while (j < rows.size() && rows[j].tx_idx == tx) {
                tx_nobs += rows[j].nobs;
                bool multi = (eq_deg[rows[j].eq_local_id] > 1);
                if (multi && !rows[j].equal) sum_noteq_rcw += rows[j].rc_width;
                ++j;
            }
            if (tx_nobs > 0.0) {
                int32_t tl = static_cast<int32_t>(active_txs.size());
                active_txs.push_back(tx);
                tx_to_local.emplace(tx, tl);
                for (size_t k = i; k < j; ++k) {
                    bool multi = (eq_deg[rows[k].eq_local_id] > 1);
                    double aval = 1.0;
                    if (d_mode) {
                        if (multi) {
                            aval = rows[k].equal
                                ? (1.0 - sum_noteq_rcw * d_rate / 1000.0)
                                : (rows[k].rc_width * d_rate / 1000.0);
                            if (d_rate == 0.0 && eq_all_not_eq[rows[k].eq_local_id]) {
                                aval = 0.01;
                            }
                            aval = std::max(0.0, std::min(1.0, aval));
                            if (rows[k].equal) {
                                aval = std::min(1.0, std::max(aval, rows[k].rc_width * d_rate / 1000.0));
                            }
                        } else {
                            aval = 1.0;
                        }
                    }
                    rows[k].aval = aval;
                    rows[k].full_aval = rows[k].equal ? aval : 0.0;
                    rows[k].unique_aval = (!multi) ? aval : 0.0;
                }
            } else {
                for (size_t k = i; k < j; ++k) rows[k].tx_idx = -1; // mark dropped
            }
            i = j;
        }

        // Collect active eq_ls (only eqs touched by >= 1 active tx, matching filterTxRc)
        std::vector<int32_t> active_eqs;
        for (const auto &r : rows) {
            if (r.tx_idx >= 0) active_eqs.push_back(r.eq_local_id);
        }
        std::sort(active_eqs.begin(), active_eqs.end());
        active_eqs.erase(std::unique(active_eqs.begin(), active_eqs.end()), active_eqs.end());
        std::vector<int32_t> eq_to_local(n_eq, -1);
        for (size_t ei = 0; ei < active_eqs.size(); ++ei) eq_to_local[active_eqs[ei]] = static_cast<int32_t>(ei);

        auto &blk = g_blocks[gi];
        blk.gene_id = active_genes[gi];
        blk.n_tx = static_cast<int32_t>(active_txs.size());
        blk.n_eq = static_cast<int32_t>(active_eqs.size());
        blk.tr_dim = static_cast<int64_t>(blk.n_tx) * static_cast<int64_t>(blk.n_eq);
        blk.gene_k = g_nobs_tot[gi];
        blk.active_txs = std::move(active_txs);
        blk.eq_nobs.resize(blk.n_eq);
        blk.eq_y.resize(blk.n_eq);
        for (size_t ei = 0; ei < active_eqs.size(); ++ei) {
            double nb = gene_eq_nobs[gi][active_eqs[ei]];
            blk.eq_nobs[ei] = nb;
            blk.eq_y[ei] = (blk.gene_k > 0.0) ? (nb / blk.gene_k) : 0.0;
        }
        for (const auto &r : rows) {
            if (r.tx_idx >= 0) {
                blk.edges.push_back(bambu_next::GeneSparseBlock::Edge{
                    tx_to_local[r.tx_idx], eq_to_local[r.eq_local_id],
                    r.aval, r.full_aval, r.unique_aval
                });
            }
        }
    }

    // Replicate exact assignGroups(): sort active genes stably by tr_dim ascending, group cumN %/% 1000 + 1
    std::vector<size_t> order_gi(g_blocks.size());
    std::iota(order_gi.begin(), order_gi.end(), 0);
    std::stable_sort(order_gi.begin(), order_gi.end(), [&](size_t a, size_t b) {
        return g_blocks[a].tr_dim < g_blocks[b].tr_dim;
    });
    std::vector<std::vector<size_t>> em_groups;
    {
        int64_t cum_n = 0;
        int64_t prev_grp = -1;
        for (size_t gi : order_gi) {
            if (g_blocks[gi].n_tx == 0) continue;
            cum_n += g_blocks[gi].tr_dim;
            int64_t grp_id = (cum_n / 1000) + 1;
            if (grp_id != prev_grp) {
                em_groups.emplace_back();
                prev_grp = grp_id;
            }
            em_groups.back().push_back(gi);
        }
    }

    // Solve sparse EM per group in parallel via OpenMP
    std::vector<double> tx_counts(n_tx, 0.0), tx_full(n_tx, 0.0), tx_uniq(n_tx, 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(n_threads) if (em_groups.size() > 64 && n_threads > 1)
#endif
    for (size_t grpi = 0; grpi < em_groups.size(); ++grpi) {
        const auto &glist = em_groups[grpi];
        int32_t tot_m = 0;
        for (size_t gi : glist) tot_m += g_blocks[gi].n_tx;

        if (tot_m == 1) {
            const auto &blk = g_blocks[glist[0]];
            double c = 0.0, f = 0.0, u = 0.0;
            for (const auto &e : blk.edges) {
                if (e.aval > 0.0) {
                    double nb = blk.eq_nobs[e.eq_l];
                    c += nb;
                    if (e.full_aval > 0.0) f += nb;
                    if (e.uniq_aval > 0.0) u += nb;
                }
            }
            int32_t tx = blk.active_txs[0];
            // Single block update is disjoint across groups -> safe without atomics if each tx belongs to 1 gene;
            // to be 100% atomic-safe even if synthetic annotations share txid across genes:
#ifdef _OPENMP
#pragma omp atomic
#endif
            tx_counts[tx] += c;
#ifdef _OPENMP
#pragma omp atomic
#endif
            tx_full[tx] += f;
#ifdef _OPENMP
#pragma omp atomic
#endif
            tx_uniq[tx] += u;
            continue;
        }

        // Precompute row sums and initialize theta = 1.0
        std::vector<std::vector<double>> thetas(glist.size()), row_sums(glist.size()), col_sums(glist.size());
        for (size_t b = 0; b < glist.size(); ++b) {
            const auto &blk = g_blocks[glist[b]];
            thetas[b].assign(blk.n_tx, 1.0);
            row_sums[b].assign(blk.n_tx, 0.0);
            col_sums[b].resize(blk.n_eq);
            for (const auto &e : blk.edges) row_sums[b][e.tx_l] += e.aval;
        }

        double delta = 1.0;
        int iter = 0;
        while (delta > conv && iter < (maxiter - 1)) {
            iter++;
            delta = 0.0;
            for (size_t b = 0; b < glist.size(); ++b) {
                const auto &blk = g_blocks[glist[b]];
                auto &th = thetas[b];
                auto &cs = col_sums[b];
                const auto &rs = row_sums[b];

                std::fill(cs.begin(), cs.end(), 0.0);
                for (const auto &e : blk.edges) {
                    cs[e.eq_l] += th[e.tx_l] * e.aval;
                }
                for (int32_t eq = 0; eq < blk.n_eq; ++eq) {
                    cs[eq] = (cs[eq] != 0.0) ? (blk.eq_y[eq] / cs[eq]) : 0.0;
                }
                std::vector<double> num(blk.n_tx, 0.0);
                for (const auto &e : blk.edges) {
                    num[e.tx_l] += (th[e.tx_l] * e.aval) * cs[e.eq_l];
                }
                for (int32_t m = 0; m < blk.n_tx; ++m) {
                    double prev_t = th[m];
                    double next_t = (rs[m] != 0.0) ? (num[m] / rs[m]) : 0.0;
                    th[m] = next_t;
                    if (next_t > minvalue) {
                        double rel = std::abs(next_t - prev_t) / next_t;
                        if (rel > delta) delta = rel;
                    }
                }
            }
        }

        // Compute final counts, fullLengthCounts, uniqueCounts
        for (size_t b = 0; b < glist.size(); ++b) {
            const auto &blk = g_blocks[glist[b]];
            const auto &th = thetas[b];
            auto &cs = col_sums[b];
            std::fill(cs.begin(), cs.end(), 0.0);
            for (const auto &e : blk.edges) cs[e.eq_l] += th[e.tx_l] * e.aval;
            for (int32_t eq = 0; eq < blk.n_eq; ++eq) {
                cs[eq] = (cs[eq] != 0.0) ? (blk.eq_nobs[eq] / cs[eq]) : 0.0; // note: K * Y == nobs
            }
            std::vector<double> c_loc(blk.n_tx, 0.0), f_loc(blk.n_tx, 0.0), u_loc(blk.n_tx, 0.0);
            for (const auto &e : blk.edges) {
                double coeff = th[e.tx_l] * cs[e.eq_l];
                c_loc[e.tx_l] += e.aval * coeff;
                f_loc[e.tx_l] += e.full_aval * coeff;
                u_loc[e.tx_l] += e.uniq_aval * coeff;
            }
            for (int32_t m = 0; m < blk.n_tx; ++m) {
                int32_t tx = blk.active_txs[m];
#ifdef _OPENMP
#pragma omp atomic
#endif
                tx_counts[tx] += c_loc[m];
#ifdef _OPENMP
#pragma omp atomic
#endif
                tx_full[tx] += f_loc[m];
#ifdef _OPENMP
#pragma omp atomic
#endif
                tx_uniq[tx] += u_loc[m];
            }
        }
    }

    double total_cnt = 0.0;
    for (double v : tx_counts) total_cnt += v;
    for (double v : incomp_by_gene) total_cnt += v;

    Rcpp::NumericVector out_counts(n_tx), out_cpm(n_tx), out_full(n_tx), out_uniq(n_tx);
    for (int32_t t = 0; t < n_tx; ++t) {
        double c = tx_counts[t];
        double cpm = (total_cnt > 0.0) ? (c / total_cnt * 1e6) : (c / total_cnt * 1e6);
        out_counts[t] = Rf_fround(c, sig_digit);
        out_cpm[t] = Rf_fround(cpm, sig_digit);
        out_full[t] = Rf_fround(tx_full[t], sig_digit);
        out_uniq[t] = tx_uniq[t];
    }

    Rcpp::CharacterVector incomp_genes(n_gene);
    Rcpp::NumericVector incomp_vals(n_gene);
    for (int32_t g = 0; g < n_gene; ++g) {
        incomp_genes[g] = ann->gene_names[g];
        incomp_vals[g] = incomp_by_gene[g];
    }

    // Optional distTable export for returnDistTable=TRUE / trackReads=TRUE
    Rcpp::List dist_export = R_NilValue;
    if (return_dist_table) {
        // Replicate modifyIncompatibleAssignment unique(DataFrame(...)) ordering cleanly
        struct DtRow {
            std::string ann_tx_id;
            int32_t txid;
            int32_t q_idx;
            double rc;
            bool comp;
            bool eq;
            int32_t dist;
            int32_t gid;
        };
        std::vector<DtRow> rows_out;
        rows_out.reserve(n_q * 2);
        for (uint8_t step = 1; step <= 3; ++step) {
            for (int32_t q = 0; q < n_q; ++q) {
                if (per_q[q].empty() || per_q[q][0].step != step) continue;
                bool ac = (rc_any_comp[q] != 0);
                for (const auto &h : per_q[q]) {
                    int32_t gid = ann->tx_geneid[h.subject_idx];
                    rows_out.push_back(DtRow{
                        ac ? ann->tx_names[h.subject_idx] : (ann->gene_names[gid] + "_unidentified"),
                        h.subject_idx + 1,
                        q,
                        rc_read_cnt[q],
                        h.compatible,
                        h.equal,
                        h.dist,
                        gid
                    });
                }
            }
        }
        // Unique preserving exact DataFrame columns (deduplicate identical tuples if uncomp replaced multiple txids?
        // Actually txid differs so only identical tuples are rare duplicate txids; if !ac, note R keeps txid, so txid distinguishes unless identical!)
        R_xlen_t m = rows_out.size();
        Rcpp::CharacterVector dt_atx(m), dt_rcid(m), dt_gid(m);
        Rcpp::IntegerVector dt_txid(m), dt_dist(m);
        Rcpp::NumericVector dt_rc(m);
        Rcpp::LogicalVector dt_comp(m), dt_eq(m);
        for (R_xlen_t i = 0; i < m; ++i) {
            const auto &r = rows_out[i];
            dt_atx[i] = r.ann_tx_id;
            dt_txid[i] = r.txid;
            dt_rcid[i] = CHAR(STRING_ELT(q_names_r, r.q_idx));
            dt_rc[i] = r.rc;
            dt_comp[i] = r.comp;
            dt_eq[i] = r.eq;
            dt_dist[i] = r.dist;
            dt_gid[i] = ann->gene_names[r.gid];
        }
        dist_export = Rcpp::List::create(
            Rcpp::Named("annotationTxId") = dt_atx,
            Rcpp::Named("txid") = dt_txid,
            Rcpp::Named("readClassId") = dt_rcid,
            Rcpp::Named("readCount") = dt_rc,
            Rcpp::Named("compatible") = dt_comp,
            Rcpp::Named("equal") = dt_eq,
            Rcpp::Named("dist") = dt_dist,
            Rcpp::Named("GENEID") = dt_gid
        );
    }

    return Rcpp::List::create(
        Rcpp::Named("counts") = out_counts,
        Rcpp::Named("CPM") = out_cpm,
        Rcpp::Named("fullLengthCounts") = out_full,
        Rcpp::Named("uniqueCounts") = out_uniq,
        Rcpp::Named("incomp_genes") = incomp_genes,
        Rcpp::Named("incomp_counts") = incomp_vals,
        Rcpp::Named("dist_table") = dist_export
    );
    END_RCPP
}
