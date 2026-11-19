#include <Rcpp.h>

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "common/flat_csr_index.hpp"
#include "stage2/incremental_min_eq.hpp"
#include "stage3/csr_splice_dist.hpp"

namespace bambu_next {
namespace {

struct IntronCoordKey {
    int32_t seqid;
    int32_t start;
    int32_t end;

    bool operator==(const IntronCoordKey &o) const noexcept {
        return seqid == o.seqid && start == o.start && end == o.end;
    }
};

struct IntronCoordKeyHash {
    size_t operator()(const IntronCoordKey &k) const noexcept {
        uint64_t a = static_cast<uint32_t>(k.seqid);
        uint64_t b = static_cast<uint32_t>(k.start);
        uint64_t c = static_cast<uint32_t>(k.end);
        uint64_t h = (b << 32) ^ c ^ (a * 0x9e3779b97f4a7c15ULL);
        h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
        h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
        return static_cast<size_t>(h ^ (h >> 31));
    }
};

struct IntronGeneEntry {
    int8_t strand;
    int32_t gene_id;
    int32_t mult;
};

struct RefIntronBucket {
    uint8_t strand_bits = 0; // bit0: '+', bit1: '-', bit2: '*'
    std::vector<IntronGeneEntry> entries;
};

static inline uint8_t strand_bit(int8_t s) noexcept {
    if (s == 1) return 1;
    if (s == -1) return 2;
    return 4;
}

static inline bool strand_compatible(int8_t a, int8_t b) noexcept {
    return a == b || a == 0 || b == 0;
}

static inline void parse_csv_ints(const char *s, std::vector<int32_t> &out) {
    out.clear();
    if (!s || *s == '\0') return;
    const char *p = s;
    while (*p) {
        while (*p == ',' || *p == ' ') ++p;
        if (!*p) break;
        int32_t val = 0;
        auto res = std::from_chars(p, p + std::strlen(p), val);
        if (res.ec == std::errc()) {
            out.push_back(val);
            p = res.ptr;
        } else {
            ++p;
        }
    }
}

static void merge_ivs(std::vector<Interval2> &v) {
    if (v.size() <= 1) return;
    std::sort(v.begin(), v.end(), [](const Interval2 &a, const Interval2 &b) {
        if (a.start != b.start) return a.start < b.start;
        return a.end < b.end;
    });
    size_t w = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        if (w == 0 || v[i].start > v[w - 1].end + 1) {
            v[w++] = v[i];
        } else if (v[i].end > v[w - 1].end) {
            v[w - 1].end = v[i].end;
        }
    }
    v.resize(w);
}

static bool cut_compatible_merged(
    const std::vector<Interval2> &q_cut_merged,
    const std::vector<Interval2> &q_gaps,
    const FlatCSRAnnotation &ann,
    int32_t sid,
    std::vector<Interval2> &s_buf
) {
    int32_t sb = ann.tx_exon_ptr[sid];
    int32_t se = ann.tx_exon_ptr[sid + 1];
    int32_t sn = se - sb;
    if (sn <= 0 || q_cut_merged.empty()) return false;
    if (q_cut_merged.front().start < ann.tx_cut_min_start[sid] ||
        q_cut_merged.back().end > ann.tx_cut_max_end[sid]) {
        return false;
    }

    const int32_t *ss = ann.cut_exon_start.data() + sb;
    const int32_t *send = ann.cut_exon_end.data() + sb;

    // Fast path when subject cut exons are already strictly separated and sorted
    bool s_clean = true;
    for (int32_t j = 1; j < sn; ++j) {
        if (ss[j] <= send[j - 1] + 1) {
            s_clean = false;
            break;
        }
    }

    if (s_clean) {
        size_t j = 0;
        for (const auto &qi : q_cut_merged) {
            while (j < static_cast<size_t>(sn) && send[j] < qi.start) ++j;
            if (j >= static_cast<size_t>(sn) || ss[j] > qi.start || qi.end > send[j]) {
                return false;
            }
        }
        size_t sj = 0;
        for (const auto &gi : q_gaps) {
            while (sj < static_cast<size_t>(sn) && send[sj] < gi.start) ++sj;
            if (sj < static_cast<size_t>(sn) && ss[sj] <= gi.end) {
                return false;
            }
        }
        return true;
    }

    s_buf.resize(sn);
    for (int32_t j = 0; j < sn; ++j) s_buf[j] = Interval2{ss[j], send[j]};
    merge_ivs(s_buf);
    size_t j = 0;
    for (const auto &qi : q_cut_merged) {
        while (j < s_buf.size() && s_buf[j].end < qi.start) ++j;
        if (j >= s_buf.size() || s_buf[j].start > qi.start || qi.end > s_buf[j].end) {
            return false;
        }
    }
    size_t sj = 0;
    for (const auto &gi : q_gaps) {
        while (sj < s_buf.size() && s_buf[sj].end < gi.start) ++sj;
        if (sj < s_buf.size() && s_buf[sj].start <= gi.end) {
            return false;
        }
    }
    return true;
}

static bool cut_equal_exact(
    const std::vector<Interval2> &q_cut_sorted,
    const FlatCSRAnnotation &ann,
    int32_t sid,
    std::vector<Interval2> &s_buf
) {
    int32_t sb = ann.tx_exon_ptr[sid];
    int32_t se = ann.tx_exon_ptr[sid + 1];
    int32_t sn = se - sb;
    if (static_cast<int32_t>(q_cut_sorted.size()) != sn) return false;
    const int32_t *ss = ann.cut_exon_start.data() + sb;
    const int32_t *send = ann.cut_exon_end.data() + sb;

    bool s_sorted = true;
    for (int32_t j = 1; j < sn; ++j) {
        if (ss[j] < ss[j - 1] || (ss[j] == ss[j - 1] && send[j] < send[j - 1])) {
            s_sorted = false;
            break;
        }
    }
    if (s_sorted) {
        for (int32_t j = 0; j < sn; ++j) {
            if (q_cut_sorted[j].start != ss[j] || q_cut_sorted[j].end != send[j]) return false;
        }
        return true;
    }
    s_buf.resize(sn);
    for (int32_t j = 0; j < sn; ++j) s_buf[j] = Interval2{ss[j], send[j]};
    std::sort(s_buf.begin(), s_buf.end(), [](const Interval2 &a, const Interval2 &b) {
        if (a.start != b.start) return a.start < b.start;
        return a.end < b.end;
    });
    for (int32_t j = 0; j < sn; ++j) {
        if (q_cut_sorted[j].start != s_buf[j].start || q_cut_sorted[j].end != s_buf[j].end) return false;
    }
    return true;
}

} // namespace
} // namespace bambu_next

extern "C" SEXP bambu_next_build_exons_from_introns_cpp(
    SEXP q_tx_seq_1b_sexp,
    SEXP q_start_sexp,
    SEXP q_end_sexp,
    SEXP q_tx_strand_code_sexp,
    SEXP q_is_vec_sexp,
    SEXP q_ie_vec_sexp
) {
    BEGIN_RCPP
    Rcpp::IntegerVector tx_seq(q_tx_seq_1b_sexp);
    Rcpp::IntegerVector tx_s(q_start_sexp);
    Rcpp::IntegerVector tx_e(q_end_sexp);
    Rcpp::IntegerVector tx_st(q_tx_strand_code_sexp);
    Rcpp::CharacterVector is_vec(q_is_vec_sexp);
    Rcpp::CharacterVector ie_vec(q_ie_vec_sexp);

    int32_t n = static_cast<int32_t>(tx_s.size());
    Rcpp::IntegerVector tx_ex_cnts(n), tx_in_cnts(n);
    std::vector<int32_t> ex_seq, ex_s, ex_e, ex_st, ex_rk, ex_erk;
    std::vector<int32_t> in_s, in_e;
    ex_seq.reserve(static_cast<size_t>(n) * 4);
    ex_s.reserve(static_cast<size_t>(n) * 4);
    ex_e.reserve(static_cast<size_t>(n) * 4);
    ex_st.reserve(static_cast<size_t>(n) * 4);
    ex_rk.reserve(static_cast<size_t>(n) * 4);
    ex_erk.reserve(static_cast<size_t>(n) * 4);
    in_s.reserve(static_cast<size_t>(n) * 3);
    in_e.reserve(static_cast<size_t>(n) * 3);

    std::vector<int32_t> istarts, iends;
    for (int32_t i = 0; i < n; ++i) {
        istarts.clear();
        iends.clear();
        SEXP is_el = STRING_ELT(is_vec, i);
        SEXP ie_el = STRING_ELT(ie_vec, i);
        if (is_el != NA_STRING && ie_el != NA_STRING) {
            bambu_next::parse_csv_ints(CHAR(is_el), istarts);
            bambu_next::parse_csv_ints(CHAR(ie_el), iends);
        }
        int32_t k = static_cast<int32_t>(std::min(istarts.size(), iends.size()));
        int32_t n_ex = k + 1;
        tx_ex_cnts[i] = n_ex;
        tx_in_cnts[i] = k;

        int32_t sid = tx_seq[i];
        int32_t st = tx_st[i];
        int32_t s0 = tx_s[i];
        int32_t e0 = tx_e[i];

        for (int32_t m = 0; m < k; ++m) {
            in_s.push_back(istarts[m]);
            in_e.push_back(iends[m]);
        }

        for (int32_t m = 0; m < n_ex; ++m) {
            int32_t cur_s = (m == 0) ? s0 : (iends[m - 1] + 1);
            int32_t cur_e = (m == k) ? e0 : (istarts[m] - 1);
            ex_seq.push_back(sid);
            ex_s.push_back(cur_s);
            ex_e.push_back(cur_e);
            ex_st.push_back(st);
            if (st != -1) {
                ex_rk.push_back(m + 1);
                ex_erk.push_back(n_ex - m);
            } else {
                ex_rk.push_back(n_ex - m);
                ex_erk.push_back(m + 1);
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("tx_ex_counts") = tx_ex_cnts,
        Rcpp::Named("ex_seq_1b") = Rcpp::wrap(ex_seq),
        Rcpp::Named("ex_start") = Rcpp::wrap(ex_s),
        Rcpp::Named("ex_end") = Rcpp::wrap(ex_e),
        Rcpp::Named("ex_strand_code") = Rcpp::wrap(ex_st),
        Rcpp::Named("ex_rank") = Rcpp::wrap(ex_rk),
        Rcpp::Named("ex_end_rank") = Rcpp::wrap(ex_erk),
        Rcpp::Named("tx_in_counts") = tx_in_cnts,
        Rcpp::Named("in_start") = Rcpp::wrap(in_s),
        Rcpp::Named("in_end") = Rcpp::wrap(in_e)
    );
    END_RCPP
}

extern "C" SEXP bambu_next_add_new_spliced_cpp(
    SEXP csr_ptr_sexp,
    SEXP q_seq_names_sexp,
    SEXP q_tx_ex_counts_sexp,
    SEXP q_ex_seq_1b_sexp,
    SEXP q_ex_start_sexp,
    SEXP q_ex_end_sexp,
    SEXP q_ex_strand_sexp,
    SEXP q_ex_rank_sexp,
    SEXP q_ex_end_rank_sexp,
    SEXP q_in_counts_sexp,
    SEXP q_in_start_sexp,
    SEXP q_in_end_sexp,
    SEXP min_exon_distance_sexp,
    SEXP min_ps_dist_sexp,
    SEXP min_ps_dist_se_sexp,
    SEXP n_threads_sexp
) {
    BEGIN_RCPP
    Rcpp::IntegerVector q_cnts_r(q_tx_ex_counts_sexp);
    int32_t n_q = static_cast<int32_t>(q_cnts_r.size());

    if (n_q == 0 || Rf_isNull(csr_ptr_sexp)) {
        Rcpp::CharacterVector out_tx(n_q, NA_STRING), out_gene(n_q, NA_STRING), out_rct(n_q, "allNew");
        Rcpp::LogicalVector out_novel(n_q, TRUE);
        return Rcpp::List::create(
            Rcpp::Named("TXNAME") = out_tx,
            Rcpp::Named("GENEID") = out_gene,
            Rcpp::Named("readClassType") = out_rct,
            Rcpp::Named("novelTranscript") = out_novel
        );
    }

    Rcpp::XPtr<bambu_next::FlatCSRAnnotation> ann(csr_ptr_sexp);
    if (ann->n_tx == 0) {
        Rcpp::CharacterVector out_tx(n_q, NA_STRING), out_gene(n_q, NA_STRING), out_rct(n_q, "allNew");
        Rcpp::LogicalVector out_novel(n_q, TRUE);
        return Rcpp::List::create(
            Rcpp::Named("TXNAME") = out_tx,
            Rcpp::Named("GENEID") = out_gene,
            Rcpp::Named("readClassType") = out_rct,
            Rcpp::Named("novelTranscript") = out_novel
        );
    }

    Rcpp::CharacterVector q_snames_r(q_seq_names_sexp);
    Rcpp::IntegerVector q_es_r(q_ex_seq_1b_sexp);
    Rcpp::IntegerVector q_s_r(q_ex_start_sexp);
    Rcpp::IntegerVector q_e_r(q_ex_end_sexp);
    Rcpp::IntegerVector q_st_r(q_ex_strand_sexp);
    Rcpp::IntegerVector q_rk_r(q_ex_rank_sexp);
    Rcpp::IntegerVector q_erk_r(q_ex_end_rank_sexp);
    Rcpp::IntegerVector q_inc_r(q_in_counts_sexp);
    Rcpp::IntegerVector q_ins_r(q_in_start_sexp);
    Rcpp::IntegerVector q_ine_r(q_in_end_sexp);

    int32_t max_dist = Rf_asInteger(min_exon_distance_sexp);
    int32_t ps_dist = Rf_asInteger(min_ps_dist_sexp);
    int32_t ps_dist_se = Rf_asInteger(min_ps_dist_se_sexp);
    int n_threads = std::max(1, Rf_asInteger(n_threads_sexp));

    std::vector<std::string> q_snames(q_snames_r.size());
    std::vector<int32_t> q_to_csr_seq(q_snames_r.size(), -1);
    for (R_xlen_t i = 0; i < q_snames_r.size(); ++i) {
        q_snames[i] = Rcpp::as<std::string>(q_snames_r[i]);
        auto it = ann->seq_to_id.find(q_snames[i]);
        if (it != ann->seq_to_id.end()) q_to_csr_seq[i] = it->second;
    }

    std::vector<int32_t> q_ptr(n_q + 1), q_in_ptr(n_q + 1);
    int32_t off = 0, in_off = 0;
    for (int32_t q = 0; q < n_q; ++q) {
        q_ptr[q] = off;
        off += q_cnts_r[q];
        q_in_ptr[q] = in_off;
        in_off += q_inc_r[q];
    }
    q_ptr[n_q] = off;
    q_in_ptr[n_q] = in_off;

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
    std::vector<int32_t> q_ins(in_off), q_ine(in_off);
    for (int32_t i = 0; i < in_off; ++i) {
        q_ins[i] = q_ins_r[i];
        q_ine[i] = q_ine_r[i];
    }

    // 1. Build reference intron hash map matching unlistedIntronsAnnotations <- unlist(myGaps(ann))
    std::unordered_map<bambu_next::IntronCoordKey, bambu_next::RefIntronBucket, bambu_next::IntronCoordKeyHash> ref_introns;
    ref_introns.reserve(static_cast<size_t>(ann->n_tx) * 4 + 64);
    std::vector<bambu_next::Interval2> ann_ex_buf;
    for (int32_t s = 0; s < ann->n_tx; ++s) {
        int32_t sb = ann->tx_exon_ptr[s];
        int32_t se = ann->tx_exon_ptr[s + 1];
        int32_t sn = se - sb;
        if (sn < 2) continue;
        int32_t csid = ann->tx_seqid[s];
        if (csid < 0) continue;
        int8_t sst = ann->tx_strand[s];
        int32_t gid = ann->tx_geneid[s];

        bool clean_sorted = true;
        for (int32_t j = sb + 1; j < se; ++j) {
            if (ann->exon_start[j] <= ann->exon_end[j - 1] + 1) {
                clean_sorted = false;
                break;
            }
        }
        auto add_intron = [&](int32_t gs, int32_t ge) {
            if (gs > ge) return;
            auto &bk = ref_introns[bambu_next::IntronCoordKey{csid, gs, ge}];
            bk.strand_bits |= bambu_next::strand_bit(sst);
            for (auto &ent : bk.entries) {
                if (ent.strand == sst && ent.gene_id == gid) {
                    ent.mult += 1;
                    return;
                }
            }
            bk.entries.push_back(bambu_next::IntronGeneEntry{sst, gid, 1});
        };
        if (clean_sorted) {
            for (int32_t j = sb + 1; j < se; ++j) {
                add_intron(ann->exon_end[j - 1] + 1, ann->exon_start[j] - 1);
            }
        } else {
            ann_ex_buf.resize(sn);
            for (int32_t j = 0; j < sn; ++j) {
                ann_ex_buf[j] = bambu_next::Interval2{ann->exon_start[sb + j], ann->exon_end[sb + j]};
            }
            bambu_next::merge_ivs(ann_ex_buf);
            for (size_t j = 1; j < ann_ex_buf.size(); ++j) {
                add_intron(ann_ex_buf[j - 1].end + 1, ann_ex_buf[j].start - 1);
            }
        }
    }

    // Per-query classification state
    std::vector<int32_t> eq_sub(n_q, -1);
    std::vector<int32_t> comp_sub(n_q, -1);
    std::vector<int32_t> intron_gene(n_q, -1);
    std::vector<uint8_t> fl_equal(n_q, 0), fl_comp(n_q, 0), fl_nwithin(n_q, 0);
    std::vector<uint8_t> fl_nlj(n_q, 0), fl_nfj(n_q, 0), fl_nj(n_q, 0), fl_allnew(n_q, 0);
    std::vector<uint8_t> fl_nfe(n_q, 0), fl_nle(n_q, 0);

    bool any_global_intron_ov = false;

    int act_thr = std::min(std::max(1, n_threads), 16);
    if (n_q < 256) act_thr = 1;

#ifdef _OPENMP
#pragma omp parallel num_threads(act_thr) reduction(||:any_global_intron_ov) if (act_thr > 1)
#endif
    {
        std::vector<int32_t> cands;
        std::vector<bambu_next::Interval2> q_cut_sorted, q_cut_merged, q_gaps, s_buf;
        std::vector<std::pair<int32_t, int32_t>> gene_cnts;
        cands.reserve(64);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
        for (int32_t q = 0; q < n_q; ++q) {
            int32_t qb = q_ptr[q], qe = q_ptr[q + 1];
            if (qb < qe) {
                int32_t csid = q_to_csr_seq[q_es[qb]];
                int8_t qst = q_st[qb];
                if (csid >= 0) {
                    int32_t qn = qe - qb;
                    q_cut_sorted.resize(qn);
                    for (int32_t k = qb; k < qe; ++k) {
                        int32_t cs = q_s[k], ce = q_e[k];
                        bool is_st = (qst != -1) ? (q_rk[k] == 1) : (q_erk[k] == 1);
                        if (is_st) cs = ce - 1;
                        bool is_en = (qst == -1) ? (q_rk[k] == 1) : (q_erk[k] == 1);
                        if (is_en) ce = cs + 1;
                        q_cut_sorted[k - qb] = bambu_next::Interval2{cs, ce};
                    }
                    std::sort(q_cut_sorted.begin(), q_cut_sorted.end(), [](const auto &a, const auto &b) {
                        if (a.start != b.start) return a.start < b.start;
                        return a.end < b.end;
                    });
                    q_cut_merged = q_cut_sorted;
                    bambu_next::merge_ivs(q_cut_merged);
                    q_gaps.clear();
                    for (size_t i = 1; i < q_cut_merged.size(); ++i) {
                        int32_t gs = q_cut_merged[i - 1].end + 1;
                        int32_t ge = q_cut_merged[i].start - 1;
                        if (gs <= ge) q_gaps.push_back(bambu_next::Interval2{gs, ge});
                    }
                    int32_t cmin = q_cut_merged.front().start;
                    int32_t cmax = q_cut_merged.back().end;

                    cands.clear();
                    ann->query_buckets(
                        ann->cut_span_index, csid, qst, false, cmin, cmax,
                        [&](int32_t tx) { cands.push_back(tx); });
                    std::sort(cands.begin(), cands.end());
                    cands.erase(std::unique(cands.begin(), cands.end()), cands.end());

                    int32_t best_eq = -1, best_cp = -1;
                    for (int32_t s : cands) {
                        if (best_eq < 0 && bambu_next::cut_equal_exact(q_cut_sorted, *ann, s, s_buf)) {
                            best_eq = s;
                        }
                        if (best_cp < 0 && bambu_next::cut_compatible_merged(q_cut_merged, q_gaps, *ann, s, s_buf)) {
                            best_cp = s;
                        }
                        if (best_eq >= 0 && best_cp >= 0) break;
                    }
                    eq_sub[q] = best_eq;
                    comp_sub[q] = best_cp;
                    if (best_eq >= 0) {
                        fl_equal[q] = 1;
                    } else if (best_cp >= 0) {
                        fl_comp[q] = 1;
                    }

                    // 2. Intron matching (updateWIntronMatches)
                    int32_t ib = q_in_ptr[q], ie = q_in_ptr[q + 1];
                    int32_t ni = ie - ib;
                    if (ni > 0) {
                        uint8_t qbit = bambu_next::strand_bit(qst);
                        bool all_m = true, any_m = false;
                        bool first_m = false, last_m = false;
                        int32_t n_unmatched = 0;
                        gene_cnts.clear();

                        for (int32_t k = ib; k < ie; ++k) {
                            auto it = ref_introns.find(bambu_next::IntronCoordKey{csid, q_ins[k], q_ine[k]});
                            bool m_exact = false;
                            if (it != ref_introns.end()) {
                                const auto &bk = it->second;
                                if ((bk.strand_bits & qbit) != 0) m_exact = true;
                                for (const auto &ent : bk.entries) {
                                    if (bambu_next::strand_compatible(qst, ent.strand)) {
                                        any_global_intron_ov = true;
                                        bool found_g = false;
                                        for (auto &gc : gene_cnts) {
                                            if (gc.first == ent.gene_id) {
                                                gc.second += ent.mult;
                                                found_g = true;
                                                break;
                                            }
                                        }
                                        if (!found_g) gene_cnts.emplace_back(ent.gene_id, ent.mult);
                                    }
                                }
                            }
                            if (m_exact) any_m = true;
                            else {
                                all_m = false;
                                n_unmatched++;
                            }
                            if (k == ib) first_m = m_exact;
                            if (k == ie - 1) last_m = m_exact;
                        }

                        if (all_m && !fl_comp[q] && !fl_equal[q]) fl_nwithin[q] = 1;
                        if (any_m && ((qst == 1 && !last_m) || (qst == -1 && !first_m))) fl_nlj[q] = 1;
                        if (any_m && ((qst == 1 && !first_m) || (qst == -1 && !last_m))) fl_nfj[q] = 1;
                        // Replicate R operator precedence in `sum(!x) > !firstJunctionMatch + !lastJunctionMatch`
                        // where `+` binds tighter than unary `!`, i.e., `!(firstJunctionMatch + !lastJunctionMatch)`.
                        int32_t unm_rhs = (!((first_m ? 1 : 0) + (!last_m ? 1 : 0))) ? 1 : 0;
                        if (any_m && n_unmatched > unm_rhs) fl_nj[q] = 1;
                        if (!any_m) fl_allnew[q] = 1;

                        if (!gene_cnts.empty()) {
                            int32_t best_g = gene_cnts[0].first;
                            int32_t max_c = gene_cnts[0].second;
                            for (size_t idx = 1; idx < gene_cnts.size(); ++idx) {
                                int32_t cg = gene_cnts[idx].first;
                                int32_t cc = gene_cnts[idx].second;
                                if (cc > max_c || (cc == max_c && ann->gene_names[cg] < ann->gene_names[best_g])) {
                                    max_c = cc;
                                    best_g = cg;
                                }
                            }
                            intron_gene[q] = best_g;
                        }
                    } else {
                        if (!fl_comp[q] && !fl_equal[q]) fl_nwithin[q] = 1;
                        fl_allnew[q] = 1;
                    }
                } else {
                    fl_allnew[q] = 1;
                }
            }
        }
    }

    // 3. Assign gene/dist updates if any intron overlapped reference annotations
    if (any_global_intron_ov) {
        std::vector<std::vector<bambu_next::DistHit>> per_q;
        bambu_next::calculate_dist_to_annotation_core(
            *ann, q_snames, q_ptr, q_es, q_s, q_e, q_st, q_rk, q_erk,
            max_dist, ps_dist, ps_dist_se, false, n_threads, per_q
        );
        for (int32_t q = 0; q < n_q; ++q) {
            const auto &hits = per_q[q];
            if (hits.empty()) continue;
            bool any_sm = false, any_em = false, any_cp = false;
            for (const auto &h : hits) {
                if (h.start_match) any_sm = true;
                if (h.end_match) any_em = true;
                if (h.compatible) any_cp = true;
            }
            if (any_cp) fl_comp[q] = 1;
            if (!any_sm && fl_nfj[q]) fl_nfe[q] = 1;
            if (!any_em && fl_nlj[q]) fl_nle[q] = 1;
        }
    }

    // 4. Pack exact R vectors matching addNewSplicedReadClasses
    Rcpp::CharacterVector out_tx(n_q), out_gene(n_q), out_rct(n_q);
    Rcpp::LogicalVector out_novel(n_q);
    std::string rct_buf;
    rct_buf.reserve(96);

    for (int32_t q = 0; q < n_q; ++q) {
        if (eq_sub[q] >= 0) {
            int32_t es = eq_sub[q];
            out_tx[q] = ann->tx_names[es];
            out_gene[q] = ann->gene_names[ann->tx_geneid[es]];
            out_novel[q] = FALSE;
        } else {
            SET_STRING_ELT(out_tx, q, NA_STRING);
            if (comp_sub[q] >= 0) {
                out_gene[q] = ann->gene_names[ann->tx_geneid[comp_sub[q]]];
            } else if (intron_gene[q] >= 0) {
                out_gene[q] = ann->gene_names[intron_gene[q]];
            } else {
                SET_STRING_ELT(out_gene, q, NA_STRING);
            }
            out_novel[q] = TRUE;
        }

        rct_buf.clear();
        auto add_tok = [&](const char *tok) {
            if (!rct_buf.empty()) rct_buf.push_back(':');
            rct_buf.append(tok);
        };
        if (fl_equal[q]) add_tok("equal");
        if (fl_comp[q]) add_tok("compatible");
        if (fl_nwithin[q]) add_tok("newWithin");
        if (fl_nlj[q]) add_tok("newLastJunction");
        if (fl_nfj[q]) add_tok("newFirstJunction");
        if (fl_nj[q]) add_tok("newJunction");
        if (fl_allnew[q]) add_tok("allNew");
        if (fl_nfe[q]) add_tok("newFirstExon");
        if (fl_nle[q]) add_tok("newLastExon");
        out_rct[q] = rct_buf;
    }

    return Rcpp::List::create(
        Rcpp::Named("TXNAME") = out_tx,
        Rcpp::Named("GENEID") = out_gene,
        Rcpp::Named("readClassType") = out_rct,
        Rcpp::Named("novelTranscript") = out_novel
    );
    END_RCPP
}

extern "C" SEXP bambu_next_filter_unspliced_overlap_cpp(
    SEXP csr_ptr_sexp,
    SEXP q_seq_names_sexp,
    SEXP q_seq_1b_sexp,
    SEXP q_start_sexp,
    SEXP q_end_sexp,
    SEXP q_strand_code_sexp,
    SEXP min_overlap_sexp
) {
    BEGIN_RCPP
    Rcpp::IntegerVector q_seq(q_seq_1b_sexp);
    Rcpp::IntegerVector q_s(q_start_sexp);
    Rcpp::IntegerVector q_e(q_end_sexp);
    Rcpp::IntegerVector q_st(q_strand_code_sexp);
    int32_t min_ov = Rf_asInteger(min_overlap_sexp);
    int32_t n_q = static_cast<int32_t>(q_s.size());

    if (n_q == 0) return Rcpp::IntegerVector(0);
    if (Rf_isNull(csr_ptr_sexp)) {
        Rcpp::IntegerVector keep(n_q);
        for (int32_t i = 0; i < n_q; ++i) keep[i] = i + 1;
        return keep;
    }
    Rcpp::XPtr<bambu_next::FlatCSRAnnotation> ann(csr_ptr_sexp);
    if (ann->n_tx == 0) {
        Rcpp::IntegerVector keep(n_q);
        for (int32_t i = 0; i < n_q; ++i) keep[i] = i + 1;
        return keep;
    }

    Rcpp::CharacterVector q_snames_r(q_seq_names_sexp);
    std::vector<int32_t> q_to_csr(q_snames_r.size(), -1);
    for (R_xlen_t i = 0; i < q_snames_r.size(); ++i) {
        auto it = ann->seq_to_id.find(Rcpp::as<std::string>(q_snames_r[i]));
        if (it != ann->seq_to_id.end()) q_to_csr[i] = it->second;
    }

    std::vector<int32_t> keep;
    keep.reserve(n_q);
    std::vector<bambu_next::Interval2> s_iv;
    for (int32_t q = 0; q < n_q; ++q) {
        int32_t csid = q_to_csr[q_seq[q] - 1];
        if (csid < 0) {
            keep.push_back(q + 1);
            continue;
        }
        int32_t qs = q_s[q], qe = q_e[q];
        int8_t st = static_cast<int8_t>(q_st[q]);
        if (qe - qs + 1 < min_ov) {
            keep.push_back(q + 1);
            continue;
        }
        bool hit = false;
        ann->query_buckets(
            ann->span_index, csid, st, false, qs, qe,
            [&](int32_t sid) {
                if (hit) return;
                int32_t sb = ann->tx_exon_ptr[sid], se = ann->tx_exon_ptr[sid + 1];
                int32_t sn = se - sb;
                int64_t ov_sum = 0;
                // Check sorted disjointness; if disjoint, sum direct exon overlaps, else merge
                bool disj = true;
                for (int32_t j = sb + 1; j < se; ++j) {
                    if (ann->exon_start[j] <= ann->exon_end[j - 1]) {
                        disj = false;
                        break;
                    }
                }
                if (disj) {
                    for (int32_t j = sb; j < se; ++j) {
                        int32_t os = std::max(qs, ann->exon_start[j]);
                        int32_t oe = std::min(qe, ann->exon_end[j]);
                        if (os <= oe) ov_sum += (oe - os + 1);
                    }
                } else {
                    s_iv.resize(sn);
                    for (int32_t j = 0; j < sn; ++j) {
                        s_iv[j] = bambu_next::Interval2{ann->exon_start[sb + j], ann->exon_end[sb + j]};
                    }
                    bambu_next::merge_ivs(s_iv);
                    for (const auto &iv : s_iv) {
                        int32_t os = std::max(qs, iv.start);
                        int32_t oe = std::min(qe, iv.end);
                        if (os <= oe) ov_sum += (oe - os + 1);
                    }
                }
                if (ov_sum >= min_ov) hit = true;
            });
        if (!hit) keep.push_back(q + 1);
    }
    return Rcpp::wrap(keep);
    END_RCPP
}

extern "C" SEXP bambu_next_extend_min_eq_csr_cpp(
    SEXP ref_csr_ptr_sexp,
    SEXP n_novel_sexp,
    SEXP seq_names_sexp,
    SEXP seq_lengths_sexp,
    SEXP tx_names_sexp,
    SEXP tx_gene_names_sexp,
    SEXP tx_seq_1b_sexp,
    SEXP tx_strand_code_sexp,
    SEXP tx_ex_counts_sexp,
    SEXP ex_start_sexp,
    SEXP ex_end_sexp,
    SEXP ex_rank_sexp,
    SEXP ex_end_rank_sexp,
    SEXP n_threads_sexp
) {
    BEGIN_RCPP
    auto *ann = new bambu_next::FlatCSRAnnotation();
    Rcpp::CharacterVector seq_nms(seq_names_sexp);
    Rcpp::IntegerVector seq_lens(seq_lengths_sexp);
    ann->n_seq = static_cast<int32_t>(seq_nms.size());
    ann->seq_names.resize(ann->n_seq);
    ann->seq_lengths.resize(ann->n_seq);
    for (int32_t i = 0; i < ann->n_seq; ++i) {
        std::string s = Rcpp::as<std::string>(seq_nms[i]);
        ann->seq_names[i] = s;
        ann->seq_to_id.emplace(s, i);
        ann->seq_lengths[i] = (i < seq_lens.size() && seq_lens[i] != NA_INTEGER) ? seq_lens[i] : -1;
    }

    Rcpp::CharacterVector tx_nms(tx_names_sexp);
    Rcpp::CharacterVector gene_nms(tx_gene_names_sexp);
    Rcpp::IntegerVector tx_seq(tx_seq_1b_sexp);
    Rcpp::IntegerVector tx_st(tx_strand_code_sexp);
    Rcpp::IntegerVector ex_cnts(tx_ex_counts_sexp);

    int32_t n_tx = static_cast<int32_t>(tx_nms.size());
    ann->n_tx = n_tx;
    ann->tx_names.resize(n_tx);
    ann->tx_seqid.resize(n_tx);
    ann->tx_strand.resize(n_tx);
    ann->tx_geneid.resize(n_tx);
    ann->tx_exon_ptr.resize(n_tx + 1);

    int32_t offset = 0;
    for (int32_t i = 0; i < n_tx; ++i) {
        ann->tx_names[i] = Rcpp::as<std::string>(tx_nms[i]);
        ann->tx_seqid[i] = tx_seq[i] - 1;
        ann->tx_strand[i] = static_cast<int8_t>(tx_st[i]);
        std::string g = Rcpp::as<std::string>(gene_nms[i]);
        auto git = ann->gene_to_id.find(g);
        if (git == ann->gene_to_id.end()) {
            int32_t gid = ann->n_gene++;
            ann->gene_to_id.emplace(g, gid);
            ann->gene_names.push_back(std::move(g));
            ann->tx_geneid[i] = gid;
        } else {
            ann->tx_geneid[i] = git->second;
        }
        ann->tx_exon_ptr[i] = offset;
        offset += ex_cnts[i];
    }
    ann->tx_exon_ptr[n_tx] = offset;

    Rcpp::IntegerVector ex_s(ex_start_sexp);
    Rcpp::IntegerVector ex_e(ex_end_sexp);
    Rcpp::IntegerVector ex_r(ex_rank_sexp);
    Rcpp::IntegerVector ex_er(ex_end_rank_sexp);
    ann->exon_start.assign(ex_s.begin(), ex_s.end());
    ann->exon_end.assign(ex_e.begin(), ex_e.end());
    ann->exon_rank.assign(ex_r.begin(), ex_r.end());
    ann->exon_end_rank.assign(ex_er.begin(), ex_er.end());

    ann->rebuild_derived_indices();

    int32_t n_novel = std::max(0, Rf_asInteger(n_novel_sexp));
    int n_threads = std::max(1, Rf_asInteger(n_threads_sexp));

    bambu_next::FlatCSRAnnotation dummy_ref;
    const bambu_next::FlatCSRAnnotation *ref_ptr = &dummy_ref;
    if (!Rf_isNull(ref_csr_ptr_sexp)) {
        Rcpp::XPtr<bambu_next::FlatCSRAnnotation> ref_xptr(ref_csr_ptr_sexp);
        if (ref_xptr->n_tx > 0 && static_cast<int32_t>(ref_xptr->eq_ptr.size()) == ref_xptr->n_tx + 1 &&
            n_novel + ref_xptr->n_tx == n_tx) {
            ref_ptr = ref_xptr.get();
        } else {
            n_novel = n_tx;
        }
    } else {
        n_novel = n_tx;
    }

    bambu_next::compute_incremental_min_eq_classes(*ann, *ref_ptr, n_novel, n_threads);

    Rcpp::XPtr<bambu_next::FlatCSRAnnotation> ptr(ann, true);
    Rcpp::IntegerVector out_vals(ann->eq_tx.size());
    Rcpp::IntegerVector out_lens(n_tx);
    for (size_t i = 0; i < ann->eq_tx.size(); ++i) out_vals[i] = ann->eq_tx[i] + 1;
    for (int32_t i = 0; i < n_tx; ++i) out_lens[i] = ann->eq_ptr[i + 1] - ann->eq_ptr[i];

    return Rcpp::List::create(
        Rcpp::Named("csr_ptr") = ptr,
        Rcpp::Named("eq_vals") = out_vals,
        Rcpp::Named("eq_lens") = out_lens
    );
    END_RCPP
}
