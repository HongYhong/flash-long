#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "common/flat_csr_index.hpp"
#include "stage1/stage1_types.hpp"

namespace bambu_next {

namespace {

struct IntronLookupKey {
    int32_t seq_id; // bsd seq_id
    int32_t start;
    int32_t end;
    int8_t strand;  // 1, -1, 0
    bool operator==(const IntronLookupKey &o) const {
        return seq_id == o.seq_id && start == o.start && end == o.end && strand == o.strand;
    }
};

struct IntronLookupHash {
    size_t operator()(const IntronLookupKey &k) const {
        uint64_t x = static_cast<uint64_t>(static_cast<uint32_t>(k.seq_id));
        x = (x * 0x9e3779b97f4a7c15ULL) ^ static_cast<uint32_t>(k.start);
        x = (x * 0xbf58476d1ce4e5b9ULL) ^ static_cast<uint32_t>(k.end);
        x = (x * 0x94d049bb133111ebULL) ^ static_cast<uint32_t>(k.strand + 2);
        return static_cast<size_t>(x ^ (x >> 30));
    }
};

struct AnnotatedIntronIndex {
    std::unordered_set<IntronLookupKey, IntronLookupHash> exact_set;

    void build(const FlatCSRAnnotation *ann, const std::vector<std::string> &bsd_seq_names) {
        exact_set.clear();
        if (!ann || ann->n_tx == 0) return;

        std::unordered_map<std::string, int32_t> bsd_seq_map;
        bsd_seq_map.reserve(bsd_seq_names.size() * 2 + 8);
        for (size_t i = 0; i < bsd_seq_names.size(); ++i) {
            bsd_seq_map.emplace(bsd_seq_names[i], static_cast<int32_t>(i));
        }
        std::vector<int32_t> ann_to_bsd(ann->n_seq, -1);
        for (int32_t i = 0; i < ann->n_seq; ++i) {
            auto it = bsd_seq_map.find(ann->seq_names[i]);
            if (it != bsd_seq_map.end()) ann_to_bsd[i] = it->second;
        }

        exact_set.reserve(ann->exon_start.size() / 2 + 64);
        for (int32_t tx = 0; tx < ann->n_tx; ++tx) {
            int32_t b = ann->tx_exon_ptr[tx];
            int32_t e = ann->tx_exon_ptr[tx + 1];
            if (e - b < 2) continue;
            int32_t sid = ann_to_bsd[ann->tx_seqid[tx]];
            if (sid < 0) continue;
            int8_t st = ann->tx_strand[tx];
            for (int32_t k = b; k + 1 < e; ++k) {
                int32_t ist = ann->exon_end[k] + 1;
                int32_t ien = ann->exon_start[k + 1] - 1;
                if (ist <= ien) {
                    exact_set.insert(IntronLookupKey{sid, ist, ien, st});
                }
            }
        }
    }

    // GenomicRanges::match(..., ignore.strand=FALSE) enforces exact strand equality ('*' matches '*' only)
    inline bool matches_strand(int32_t seq_id, int32_t start, int32_t end, int8_t st) const {
        return exact_set.find(IntronLookupKey{seq_id, start, end, st}) != exact_set.end();
    }
};

struct SpliceFeatureBundle {
    // Subset indices into unique_junctions that map to each row in the feature matrix (length n_all)
    // Order matching unique(GRanges(..., pos, pos)) for Start, and sort(unique(GRanges(..., pos, pos))) for End
    std::vector<int32_t> junc_to_site_idx; // for each junction j -> site row 0..n_site-1
    std::vector<int32_t> my_set_all;       // 0-based site indices matching mySet.all
    std::vector<int32_t> my_set_training;  // 0-based indices into my_set_all
    std::vector<uint8_t> train_labels;     // length my_set_all.size(), bool annotatedSplice
    std::vector<double> col_ratio;         // A.1
    std::vector<double> col_score;         // A.2
    std::vector<double> col_dist;          // A.3
    std::vector<double> col_st_prime_pos;  // A.4
    std::vector<double> col_st_prime_neg;  // A.5
    std::vector<double> col_st_pos;        // A.6 (not in model.matrix formula but preserved for completeness)
};

// Replicate createSpliceMetadata + testSpliceSites matrix generation on a subset mask of unique_junctions
static SpliceFeatureBundle build_splice_feature_bundle(
    const std::vector<UniqueJunctionRec> &uj,
    const std::vector<int32_t> &subset_indices, // sorted by (seq_id, start, end)
    bool is_start_splice,                       // TRUE="Start", FALSE="End"
    bool is_start_prime                         // TRUE="start", FALSE="end"
) {
    SpliceFeatureBundle out;
    size_t n_sub = subset_indices.size();
    out.junc_to_site_idx.assign(n_sub, -1);
    if (n_sub == 0) return out;

    // Deduplicate sites on (seq_id, pos) keeping FIRST occurrence in subset_indices (matches GRanges unique())
    struct SiteRec {
        int32_t seq_id;
        int32_t pos;
        int32_t score;
        int8_t splice_strand;
        int8_t annotated;
        int32_t orig_order;
    };
    std::vector<SiteRec> sites;
    sites.reserve(n_sub);

    auto pack_pos = [](int32_t seq_id, int32_t pos) -> uint64_t {
        return (static_cast<uint64_t>(static_cast<uint32_t>(seq_id)) << 32) |
               static_cast<uint32_t>(pos);
    };
    std::unordered_map<uint64_t, int32_t> seen;
    seen.reserve(n_sub * 2 + 16);

    for (size_t i = 0; i < n_sub; ++i) {
        const auto &rec = uj[subset_indices[i]];
        int32_t pos = is_start_splice ? rec.start : rec.end;
        uint64_t key = pack_pos(rec.seq_id, pos);
        auto it = seen.find(key);
        if (it == seen.end()) {
            int32_t sidx = static_cast<int32_t>(sites.size());
            seen.emplace(key, sidx);
            sites.push_back(SiteRec{
                rec.seq_id, pos,
                is_start_splice ? rec.start_score : rec.end_score,
                rec.splice_strand,
                is_start_splice ? rec.annotated_start : rec.annotated_end,
                sidx
            });
        }
    }

    int32_t n_site = static_cast<int32_t>(sites.size());
    if (!is_start_splice) {
        // R line 169: annotatedJunctionsTmp <- sort(unique(annotatedJunctionsTmp)) for End
        std::sort(sites.begin(), sites.end(), [](const SiteRec &a, const SiteRec &b) {
            if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
            return a.pos < b.pos;
        });
        seen.clear();
        for (int32_t s = 0; s < n_site; ++s) {
            seen.emplace(pack_pos(sites[s].seq_id, sites[s].pos), s);
        }
    }
    for (size_t i = 0; i < n_sub; ++i) {
        const auto &rec = uj[subset_indices[i]];
        int32_t pos = is_start_splice ? rec.start : rec.end;
        out.junc_to_site_idx[i] = seen[pack_pos(rec.seq_id, pos)];
    }

    for (int32_t i = 0; i < n_site; ++i) {
        const auto &cur = sites[i];
        int32_t dist_prime = 0;
        int8_t annot_prime = 0;
        int32_t score_prime = 0;
        int8_t strand_prime = 0; // 0 represents FALSE/unmatched or '*'

        if (is_start_prime) {
            if (i > 0 && sites[i - 1].seq_id == cur.seq_id) {
                dist_prime = cur.pos - sites[i - 1].pos;
                annot_prime = sites[i - 1].annotated;
                score_prime = sites[i - 1].score;
                strand_prime = sites[i - 1].splice_strand;
            }
        } else {
            if (i + 1 < n_site && sites[i + 1].seq_id == cur.seq_id) {
                dist_prime = -(sites[i + 1].pos - cur.pos);
                annot_prime = sites[i + 1].annotated;
                score_prime = sites[i + 1].score;
                strand_prime = sites[i + 1].splice_strand;
            }
        }

        bool cond_dist = is_start_prime ? (dist_prime < 15) : (std::abs(dist_prime) < 15);
        if (dist_prime != 0 && cur.splice_strand != 0 && cur.score > 0 && cond_dist) {
            int32_t m_idx = static_cast<int32_t>(out.my_set_all.size());
            out.my_set_all.push_back(i);
            if (annot_prime || cur.annotated) {
                out.my_set_training.push_back(m_idx);
            }
            out.train_labels.push_back(static_cast<uint8_t>(cur.annotated != 0 ? 1 : 0));
            double sc = static_cast<double>(cur.score);
            double scp = static_cast<double>(score_prime);
            out.col_ratio.push_back(sc / (sc + scp));
            out.col_score.push_back(sc);
            out.col_dist.push_back(static_cast<double>(dist_prime));
            out.col_st_prime_pos.push_back((strand_prime == 1) ? 1.0 : 0.0);
            out.col_st_prime_neg.push_back((strand_prime == -1) ? 1.0 : 0.0);
            out.col_st_pos.push_back((cur.splice_strand == 1) ? 1.0 : 0.0);
        }
    }
    return out;
}

static Rcpp::List bundle_to_r_matrix(const SpliceFeatureBundle &b, bool include_training) {
    size_t n = b.my_set_all.size();
    Rcpp::NumericMatrix mat(static_cast<int>(n), 6);
    for (size_t i = 0; i < n; ++i) {
        mat(i, 0) = 1.0; // (Intercept)
        mat(i, 1) = b.col_ratio[i];
        mat(i, 2) = b.col_score[i];
        mat(i, 3) = b.col_dist[i];
        mat(i, 4) = b.col_st_prime_pos[i];
        mat(i, 5) = b.col_st_prime_neg[i];
    }
    Rcpp::CharacterVector cnames = Rcpp::CharacterVector::create(
        "(Intercept)", "A.1", "A.2", "A.3", "A.4TRUE", "A.5TRUE");
    Rcpp::List dn = Rcpp::List::create(R_NilValue, cnames);
    mat.attr("dimnames") = dn;

    Rcpp::IntegerVector site_idx_1b(b.junc_to_site_idx.size());
    for (size_t i = 0; i < b.junc_to_site_idx.size(); ++i) site_idx_1b[i] = b.junc_to_site_idx[i] + 1;

    Rcpp::IntegerVector set_all_1b(n);
    for (size_t i = 0; i < n; ++i) set_all_1b[i] = b.my_set_all[i] + 1;

    if (!include_training) {
        return Rcpp::List::create(
            Rcpp::Named("mat") = mat,
            Rcpp::Named("my_set_all") = set_all_1b,
            Rcpp::Named("junc_to_site") = site_idx_1b
        );
    }

    Rcpp::IntegerVector tr_idx_1b(b.my_set_training.size());
    Rcpp::LogicalVector tr_lbl(b.my_set_training.size());
    for (size_t i = 0; i < b.my_set_training.size(); ++i) {
        int32_t m = b.my_set_training[i];
        tr_idx_1b[i] = m + 1;
        tr_lbl[i] = (b.train_labels[m] != 0);
    }
    return Rcpp::List::create(
        Rcpp::Named("mat") = mat,
        Rcpp::Named("my_set_all") = set_all_1b,
        Rcpp::Named("junc_to_site") = site_idx_1b,
        Rcpp::Named("train_idx") = tr_idx_1b,
        Rcpp::Named("train_labels") = tr_lbl
    );
}

// Replicate findJunctionsByStrand resolution on candidate subset indices (sorted by seq_id, start, end)
static void resolve_strand_junctions(
    const std::vector<UniqueJunctionRec> &uj,
    const std::vector<int32_t> &cand_indices,
    const std::vector<int32_t> &hc_subset_pos, // indices into cand_indices that were predicted
    const std::vector<double> &p_ss,
    const std::vector<double> &p_se,
    const std::vector<double> &p_es,
    const std::vector<double> &p_ee,
    std::vector<int32_t> &out_merged_hc // global index into uj, or -1 if NA
) {
    size_t n_cand = cand_indices.size();
    out_merged_hc.assign(n_cand, -1);

    std::vector<uint8_t> is_hc_ref(n_cand, 0);
    for (size_t k = 0; k < hc_subset_pos.size(); ++k) {
        int32_t cpos = hc_subset_pos[k];
        int32_t gi = cand_indices[cpos];
        if (uj[gi].annotated_junction) {
            is_hc_ref[cpos] = 1;
            continue;
        }
        int cnt = 0;
        if (std::isnan(p_ss[k]) || p_ss[k] > 0.5) cnt++;
        if (std::isnan(p_se[k]) || p_se[k] > 0.5) cnt++;
        if (std::isnan(p_es[k]) || p_es[k] > 0.5) cnt++;
        if (std::isnan(p_ee[k]) || p_ee[k] > 0.5) cnt++;
        if (cnt == 4) is_hc_ref[cpos] = 1;
    }

    // Build subject list sorted by (seq_id, start, cand_pos) and count HC subjects per seq_id
    struct SubEntry {
        int32_t seq_id;
        int32_t start;
        int32_t end;
        int32_t cand_pos;
        int32_t global_idx;
    };
    std::vector<SubEntry> subs;
    subs.reserve(hc_subset_pos.size());
    int32_t max_seq = 0;
    for (size_t c = 0; c < n_cand; ++c) {
        int32_t s = uj[cand_indices[c]].seq_id;
        if (s > max_seq) max_seq = s;
        if (is_hc_ref[c]) {
            int32_t gi = cand_indices[c];
            subs.push_back(SubEntry{uj[gi].seq_id, uj[gi].start, uj[gi].end, static_cast<int32_t>(c), gi});
        }
    }
    if (subs.empty()) return;

    std::vector<int32_t> hc_per_seq(max_seq + 1, 0);
    for (const auto &se : subs) {
        if (se.seq_id >= 0) hc_per_seq[se.seq_id]++;
    }
    std::vector<int32_t> rem_per_seq(max_seq + 1, 0);

    // Note: subs is already sorted by (seq_id, start, end, cand_pos) because cand_indices is sorted.
    // Match IRanges::findOverlaps_GNCList per-chromosome NCList traversal:
    //   n_q_chr >= n_s_chr -> NCList on subject -> traverses (start ASC, end DESC) -> last wins (max start, min end)
    //   n_q_chr <  n_s_chr -> NCList on query   -> traverses (start ASC, end ASC)  -> last wins (max start, max end)
    for (int maxDist = 0; maxDist <= 10; ++maxDist) {
        std::fill(rem_per_seq.begin(), rem_per_seq.end(), 0);
        for (size_t c = 0; c < n_cand; ++c) {
            if (out_merged_hc[c] < 0) {
                int32_t s = uj[cand_indices[c]].seq_id;
                if (s >= 0) rem_per_seq[s]++;
            }
        }
        for (size_t c = 0; c < n_cand; ++c) {
            if (out_merged_hc[c] >= 0) continue;
            const auto &q = uj[cand_indices[c]];
            bool subj_nclist = (q.seq_id >= 0 && rem_per_seq[q.seq_id] >= hc_per_seq[q.seq_id]);
            SubEntry key{q.seq_id, q.start - maxDist, 0, 0, 0};
            auto it = std::lower_bound(subs.begin(), subs.end(), key, [](const SubEntry &a, const SubEntry &b) {
                if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
                return a.start < b.start;
            });
            int32_t best_gi = -1;
            int32_t best_s = 0, best_e = 0, best_cpos = -1;
            int32_t max_sub_end = q.start + maxDist;
            for (; it != subs.end() && it->seq_id == q.seq_id && it->start <= max_sub_end; ++it) {
                if (std::abs(q.end - it->end) <= maxDist) {
                    bool better = (best_gi < 0);
                    if (!better) {
                        if (it->start != best_s) {
                            better = (it->start > best_s);
                        } else if (it->end != best_e) {
                            better = subj_nclist ? (it->end < best_e) : (it->end > best_e);
                        } else {
                            better = (it->cand_pos > best_cpos);
                        }
                    }
                    if (better) {
                        best_gi = it->global_idx;
                        best_s = it->start;
                        best_e = it->end;
                        best_cpos = it->cand_pos;
                    }
                }
            }
            if (best_gi >= 0) {
                out_merged_hc[c] = best_gi;
            }
        }
    }
}

} // namespace

// Call 2: Annotate junctions + run strand correction + prepare XGBoost training/scoring matrices
extern "C" SEXP bambu_next_prepare_junc_ec_cpp(
    SEXP bsd_ptr_sexp,
    SEXP csr_ptr_sexp,
    SEXP stranded_sexp
) {
    BEGIN_RCPP
    Rcpp::XPtr<bambu_next::BamSampleData> bsd(bsd_ptr_sexp);
    const FlatCSRAnnotation *ann = Rf_isNull(csr_ptr_sexp)
        ? nullptr
        : Rcpp::XPtr<FlatCSRAnnotation>(csr_ptr_sexp).get();
    bool stranded = Rcpp::as<bool>(stranded_sexp);

    size_t n_j = bsd->unique_junctions.size();
    if (n_j == 0) {
        return Rcpp::List::create(Rcpp::Named("empty") = true);
    }

    AnnotatedIntronIndex ann_idx;
    ann_idx.build(ann, bsd->seq_names);

    // 1. Initialize corrected_strand = splice_strand & run junctionStrandCorrection
    for (size_t j = 0; j < n_j; ++j) {
        bsd->unique_junctions[j].corrected_strand = bsd->unique_junctions[j].splice_strand;
    }

    auto count_annot_overlap = [&](const std::vector<int8_t> &st_vec) -> int32_t {
        int32_t cnt = 0;
        for (size_t j = 0; j < n_j; ++j) {
            const auto &uj = bsd->unique_junctions[j];
            if (ann_idx.matches_strand(uj.seq_id, uj.start, uj.end, st_vec[j])) cnt++;
        }
        return cnt;
    };

    std::vector<int8_t> cur_st(n_j), upd_st(n_j);
    for (size_t j = 0; j < n_j; ++j) cur_st[j] = bsd->unique_junctions[j].splice_strand;
    int32_t cur_annot_num = count_annot_overlap(cur_st);

    if (cur_annot_num > 0) {
        std::vector<int32_t> read_plus(n_j), read_minus(n_j);
        while (true) {
            if (!stranded) {
                std::fill(read_plus.begin(), read_plus.end(), 0);
                std::fill(read_minus.begin(), read_minus.end(), 0);
                int32_t n_sp = bsd->spliced.n_reads;
                for (int32_t r = 0; r < n_sp; ++r) {
                    int32_t b = bsd->spliced.junc_ptr[r];
                    int32_t e = bsd->spliced.junc_ptr[r + 1];
                    int32_t p_cnt = 0, m_cnt = 0;
                    for (int32_t k = b; k < e; ++k) {
                        int8_t js = cur_st[bsd->spliced.junc_ids[k]];
                        if (js == 1) p_cnt++;
                        else if (js == -1) m_cnt++;
                    }
                    int32_t diff = m_cnt - p_cnt;
                    if (diff < 0) {
                        for (int32_t k = b; k < e; ++k) read_plus[bsd->spliced.junc_ids[k]]++;
                    } else if (diff > 0) {
                        for (int32_t k = b; k < e; ++k) read_minus[bsd->spliced.junc_ids[k]]++;
                    }
                }
                for (size_t j = 0; j < n_j; ++j) {
                    int32_t d = read_minus[j] - read_plus[j];
                    upd_st[j] = (d < 0) ? 1 : ((d > 0) ? -1 : cur_st[j]);
                }
            } else {
                for (size_t j = 0; j < n_j; ++j) {
                    int32_t d = bsd->unique_junctions[j].minus_score - bsd->unique_junctions[j].plus_score;
                    upd_st[j] = (d < 0) ? 1 : ((d > 0) ? -1 : cur_st[j]);
                }
            }
            int32_t new_annot = count_annot_overlap(upd_st);
            if (new_annot > cur_annot_num) {
                cur_annot_num = new_annot;
                cur_st = upd_st;
            } else {
                break;
            }
        }
    }

    for (size_t j = 0; j < n_j; ++j) {
        bsd->unique_junctions[j].corrected_strand = cur_st[j];
        const auto &uj = bsd->unique_junctions[j];
        bsd->unique_junctions[j].annotated_junction =
            ann_idx.matches_strand(uj.seq_id, uj.start, uj.end, cur_st[j]) ? 1 : 0;
    }

    // 2. Per-chromosome annotatedStart and annotatedEnd (%in% start/end[annotatedJunction])
    int32_t n_annot_j = 0;
    int32_t n_unannot_j = 0;
    {
        size_t i = 0;
        std::unordered_set<int32_t> ast, aen;
        while (i < n_j) {
            size_t s = i;
            int32_t seq = bsd->unique_junctions[i].seq_id;
            ast.clear();
            aen.clear();
            while (i < n_j && bsd->unique_junctions[i].seq_id == seq) {
                if (bsd->unique_junctions[i].annotated_junction) {
                    ast.insert(bsd->unique_junctions[i].start);
                    aen.insert(bsd->unique_junctions[i].end);
                    n_annot_j++;
                } else {
                    n_unannot_j++;
                }
                i++;
            }
            for (size_t k = s; k < i; ++k) {
                bsd->unique_junctions[k].annotated_start =
                    (ast.find(bsd->unique_junctions[k].start) != ast.end()) ? 1 : 0;
                bsd->unique_junctions[k].annotated_end =
                    (aen.find(bsd->unique_junctions[k].end) != aen.end()) ? 1 : 0;
            }
        }
    }

    bool need_fit = (n_annot_j > 5000 && n_unannot_j > 4000);
    std::vector<int32_t> all_idx(n_j);
    for (size_t i = 0; i < n_j; ++i) all_idx[i] = static_cast<int32_t>(i);

    Rcpp::List train_bundles = R_NilValue;
    if (need_fit) {
        train_bundles = Rcpp::List::create(
            Rcpp::Named("Start.start") = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, all_idx, true, true), true),
            Rcpp::Named("Start.end")   = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, all_idx, true, false), true),
            Rcpp::Named("End.start")   = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, all_idx, false, true), true),
            Rcpp::Named("End.end")     = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, all_idx, false, false), true)
        );
    }

    // 3. Also prepare feature bundles for findHighConfidenceJunctions (Plus and Minus subsets)
    std::vector<int32_t> cand_plus, cand_minus, hc_plus_sub, hc_minus_sub;
    cand_plus.reserve(n_j / 2 + 16);
    cand_minus.reserve(n_j / 2 + 16);
    for (size_t j = 0; j < n_j; ++j) {
        int8_t st = bsd->unique_junctions[j].corrected_strand;
        if (st == 1 || st == 0) cand_plus.push_back(static_cast<int32_t>(j));
        if (st == -1 || st == 0) cand_minus.push_back(static_cast<int32_t>(j));
    }
    for (int32_t gi : cand_plus) {
        const auto &uj = bsd->unique_junctions[gi];
        if (uj.score > 1 && uj.splice_strand == 1) hc_plus_sub.push_back(gi);
    }
    for (int32_t gi : cand_minus) {
        const auto &uj = bsd->unique_junctions[gi];
        if (uj.score > 1 && uj.splice_strand == -1) hc_minus_sub.push_back(gi);
    }

    bool can_run_hc = (!hc_plus_sub.empty() && !hc_minus_sub.empty());
    Rcpp::List hc_plus_bundles = R_NilValue;
    Rcpp::List hc_minus_bundles = R_NilValue;
    if (can_run_hc) {
        hc_plus_bundles = Rcpp::List::create(
            Rcpp::Named("Start.start") = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, hc_plus_sub, true, true), false),
            Rcpp::Named("Start.end")   = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, hc_plus_sub, true, false), false),
            Rcpp::Named("End.start")   = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, hc_plus_sub, false, true), false),
            Rcpp::Named("End.end")     = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, hc_plus_sub, false, false), false)
        );
        hc_minus_bundles = Rcpp::List::create(
            Rcpp::Named("Start.start") = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, hc_minus_sub, true, true), false),
            Rcpp::Named("Start.end")   = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, hc_minus_sub, true, false), false),
            Rcpp::Named("End.start")   = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, hc_minus_sub, false, true), false),
            Rcpp::Named("End.end")     = bundle_to_r_matrix(build_splice_feature_bundle(bsd->unique_junctions, hc_minus_sub, false, false), false)
        );
    }

    return Rcpp::List::create(
        Rcpp::Named("empty") = false,
        Rcpp::Named("need_fit") = need_fit,
        Rcpp::Named("can_run_hc") = can_run_hc,
        Rcpp::Named("train_bundles") = train_bundles,
        Rcpp::Named("hc_plus_bundles") = hc_plus_bundles,
        Rcpp::Named("hc_minus_bundles") = hc_minus_bundles
    );
    END_RCPP
}

// Helper to map XGBoost predictions back into BamSampleData unique_junctions
void apply_junc_ec_predictions(
    BamSampleData &bsd,
    bool can_run_hc,
    const std::vector<double> &plus_ss, const std::vector<double> &plus_se,
    const std::vector<double> &plus_es, const std::vector<double> &plus_ee,
    const std::vector<double> &minus_ss, const std::vector<double> &minus_se,
    const std::vector<double> &minus_es, const std::vector<double> &minus_ee
) {
    size_t n_j = bsd.unique_junctions.size();
    if (!can_run_hc) {
        for (size_t j = 0; j < n_j; ++j) {
            bsd.unique_junctions[j].merged_hc_id = static_cast<int32_t>(j);
            bsd.unique_junctions[j].merged_hc_id_no_na = static_cast<int32_t>(j);
        }
        return;
    }

    std::vector<int32_t> cand_plus, cand_minus, hc_plus_pos, hc_minus_pos;
    std::vector<int32_t> global_to_plus_pos(n_j, -1), global_to_minus_pos(n_j, -1);
    for (size_t j = 0; j < n_j; ++j) {
        int8_t st = bsd.unique_junctions[j].corrected_strand;
        if (st == 1 || st == 0) {
            int32_t p = static_cast<int32_t>(cand_plus.size());
            global_to_plus_pos[j] = p;
            cand_plus.push_back(static_cast<int32_t>(j));
            if (bsd.unique_junctions[j].score > 1 && bsd.unique_junctions[j].splice_strand == 1) {
                hc_plus_pos.push_back(p);
            }
        }
        if (st == -1 || st == 0) {
            int32_t p = static_cast<int32_t>(cand_minus.size());
            global_to_minus_pos[j] = p;
            cand_minus.push_back(static_cast<int32_t>(j));
            if (bsd.unique_junctions[j].score > 1 && bsd.unique_junctions[j].splice_strand == -1) {
                hc_minus_pos.push_back(p);
            }
        }
    }

    std::vector<int32_t> merged_plus, merged_minus;
    resolve_strand_junctions(bsd.unique_junctions, cand_plus, hc_plus_pos,
                             plus_ss, plus_se, plus_es, plus_ee, merged_plus);
    resolve_strand_junctions(bsd.unique_junctions, cand_minus, hc_minus_pos,
                             minus_ss, minus_se, minus_es, minus_ee, merged_minus);

    // useRefJunctionForConflict + merge
    for (size_t j = 0; j < n_j; ++j) {
        int32_t pp = global_to_plus_pos[j];
        int32_t mp = global_to_minus_pos[j];
        int32_t res_p = (pp >= 0) ? merged_plus[pp] : -1;
        int32_t res_m = (mp >= 0) ? merged_minus[mp] : -1;
        int32_t final_id = -1;
        if (res_p >= 0 && res_m >= 0) {
            int32_t diff = bsd.unique_junctions[res_p].score - bsd.unique_junctions[res_m].score;
            if (diff > 0) final_id = res_p;
            else if (diff < 0) final_id = res_m;
            else final_id = -1; // conflict tie sets NA
        } else if (mp >= 0) {
            // Replicate Bambu useRefJunctionForConflict line 269 overwriting '*' entries with candidateJunctionsMinus
            final_id = res_m;
        } else {
            final_id = res_p;
        }
        bsd.unique_junctions[j].merged_hc_id = final_id;
        bsd.unique_junctions[j].merged_hc_id_no_na = (final_id >= 0) ? final_id : static_cast<int32_t>(j);
    }
}

} // namespace bambu_next
