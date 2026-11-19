#include <Rcpp.h>

#include <algorithm>
#include <cstdint>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "stage3/csr_splice_dist.hpp"

namespace bambu_next {
namespace {

static inline bool ranges_equal_sorted(const std::vector<Interval2> &a,
                                       const int32_t *sb, const int32_t *eb, int32_t sn) {
    if (static_cast<int32_t>(a.size()) != sn) return false;
    for (int32_t i = 0; i < sn; ++i) {
        if (a[i].start != sb[i] || a[i].end != eb[i]) return false;
    }
    return true;
}

static void merge_disjoint(std::vector<Interval2> &v) {
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

// Compute width(setdiff(a, b)) total sum and count of contiguous pieces >= min_dist
static void setdiff_stats(
    const std::vector<Interval2> &a,
    const std::vector<Interval2> &b,
    int32_t min_dist,
    int32_t &out_total,
    int32_t &out_ge_count
) {
    out_total = 0;
    out_ge_count = 0;
    size_t j = 0;
    for (const auto &iv : a) {
        int32_t cur = iv.start;
        int32_t end = iv.end;
        while (j < b.size() && b[j].end < cur) ++j;
        size_t jj = j;
        while (jj < b.size() && b[jj].start <= end) {
            if (b[jj].start > cur) {
                int32_t len = b[jj].start - cur;
                out_total += len;
                if (len >= min_dist) out_ge_count++;
            }
            if (b[jj].end + 1 > cur) cur = b[jj].end + 1;
            if (cur > end) break;
            ++jj;
        }
        if (cur <= end) {
            int32_t len = end - cur + 1;
            out_total += len;
            if (len >= min_dist) out_ge_count++;
        }
    }
}

// Compute width(intersect(a, b)) total sum and count of contiguous pieces >= min_dist
static void intersect_stats(
    const std::vector<Interval2> &a,
    const std::vector<Interval2> &b,
    int32_t min_dist,
    int32_t &out_total,
    int32_t &out_ge_count
) {
    out_total = 0;
    out_ge_count = 0;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        int32_t s = std::max(a[i].start, b[j].start);
        int32_t e = std::min(a[i].end, b[j].end);
        if (s <= e) {
            int32_t len = e - s + 1;
            out_total += len;
            if (len >= min_dist) out_ge_count++;
        }
        if (a[i].end < b[j].end) ++i;
        else ++j;
    }
}

// Replicate calculateFirstLastExonsDist for single query exon q_ex against subjectExon s_ex and full subject exons
static void first_last_exon_dist(
    int32_t q_s, int32_t q_e,
    int32_t s_term_s, int32_t s_term_e,
    const int32_t *subj_s, const int32_t *subj_e, int32_t sn,
    bool &out_match,
    int32_t &out_uq,
    int32_t &out_us
) {
    out_match = (q_s <= s_term_e && s_term_s <= q_e);

    // ExonDist = sum_{j}(width(subj[j])*(ov_j>0) - 2*ov_j) + width(q_ex)
    int32_t q_w = q_e - q_s + 1;
    int64_t ex_dist = q_w;
    std::vector<Interval2> cov;
    cov.reserve(4);
    for (int32_t j = 0; j < sn; ++j) {
        int32_t os = std::max(q_s, subj_s[j]);
        int32_t oe = std::min(q_e, subj_e[j]);
        if (os <= oe) {
            int32_t ov = oe - os + 1;
            int32_t sw = subj_e[j] - subj_s[j] + 1;
            ex_dist += static_cast<int64_t>(sw) - 2LL * ov;
            cov.push_back(Interval2{os, oe});
        }
    }
    merge_disjoint(cov);
    int32_t covered_q = 0;
    for (const auto &c : cov) covered_q += (c.end - c.start + 1);
    out_uq = q_w - covered_q;
    out_us = static_cast<int32_t>(ex_dist - out_uq);
}

static void score_single_hit(
    const FlatCSRAnnotation &ann,
    int32_t sid,
    int8_t q_st,
    const std::vector<Interval2> &q_cut_merged,
    const std::vector<Interval2> &q_cut_raw,
    const std::vector<Interval2> &q_splice,
    int32_t q_start_s, int32_t q_start_e,
    int32_t q_end_s, int32_t q_end_e,
    int32_t step_max_dist,
    bool dist_calculated,
    DistHit &out_hit
) {
    int32_t sb = ann.tx_exon_ptr[sid];
    int32_t se = ann.tx_exon_ptr[sid + 1];
    int32_t sn = se - sb;
    const int32_t *ss = ann.exon_start.data() + sb;
    const int32_t *send = ann.exon_end.data() + sb;

    // subject_full merged for rangesDist
    std::vector<Interval2> s_full(sn);
    for (int32_t i = 0; i < sn; ++i) s_full[i] = Interval2{ss[i], send[i]};
    merge_disjoint(s_full);

    int32_t uq = 0, q_out = 0, us = 0, s_out = 0;
    setdiff_stats(q_cut_merged, s_full, step_max_dist, uq, q_out);
    intersect_stats(s_full, q_splice, step_max_dist, us, s_out);

    // Check equal against cutStartEnd(subject)
    bool is_eq = false;
    if (q_st == ann.tx_strand[sid] || q_st == 0 || ann.tx_strand[sid] == 0) {
        is_eq = ranges_equal_sorted(q_cut_raw, ann.cut_exon_start.data() + sb,
                                    ann.cut_exon_end.data() + sb, sn);
    }

    // Subject start & end exon (exon_rank == 1 and exon_endRank == 1)
    int32_t s_st_s = ss[0], s_st_e = send[0];
    int32_t s_en_s = ss[0], s_en_e = send[0];
    for (int32_t k = sb; k < se; ++k) {
        if (ann.exon_rank[k] == 1) {
            s_st_s = ann.exon_start[k];
            s_st_e = ann.exon_end[k];
        }
        if (ann.exon_end_rank[k] == 1) {
            s_en_s = ann.exon_start[k];
            s_en_e = ann.exon_end[k];
        }
    }

    bool sm = false, em = false;
    int32_t usq = 0, uss = 0, ueq = 0, ues = 0;
    first_last_exon_dist(q_start_s, q_start_e, s_st_s, s_st_e, ss, send, sn, sm, usq, uss);
    first_last_exon_dist(q_end_s, q_end_e, s_en_s, s_en_e, ss, send, sn, em, ueq, ues);

    out_hit.subject_idx = sid;
    out_hit.uq = uq;
    out_hit.us = us;
    out_hit.q_out = q_out;
    out_hit.s_out = s_out;
    out_hit.compatible = (q_out == 0 && s_out == 0);
    out_hit.equal = is_eq;
    out_hit.start_match = sm;
    out_hit.end_match = em;
    out_hit.u_start_q = usq;
    out_hit.u_start_s = uss;
    out_hit.u_end_q = ueq;
    out_hit.u_end_s = ues;
    out_hit.dist = dist_calculated ? (uq + us + usq + ueq) : (uq + us);
    out_hit.step = 1;
}

static void filter_hits_for_query(
    std::vector<DistHit> &hits,
    int32_t ps_dist,
    int32_t ps_se,
    bool dist_calculated,
    bool use_nclist_subj
) {
    if (hits.empty()) return;
    int32_t min_d = hits[0].dist;
    for (const auto &h : hits) if (h.dist < min_d) min_d = h.dist;

    std::vector<DistHit> step1;
    step1.reserve(hits.size());
    for (const auto &h : hits) {
        if (h.dist <= min_d + ps_dist) step1.push_back(h);
    }

    if (!dist_calculated && !step1.empty()) {
        int32_t min_out = step1[0].q_out + step1[0].s_out;
        for (const auto &h : step1) {
            int32_t s = h.q_out + h.s_out;
            if (s < min_out) min_out = s;
        }
        std::vector<DistHit> step2;
        step2.reserve(step1.size());
        for (const auto &h : step1) {
            if (h.q_out + h.s_out == min_out) step2.push_back(h);
        }
        int best_se = 0;
        for (const auto &h : step2) {
            if (h.u_start_q <= ps_se && h.u_end_q <= ps_se) {
                best_se = 1;
                break;
            }
        }
        step1.clear();
        for (const auto &h : step2) {
            int cur_se = (h.u_start_q <= ps_se && h.u_end_q <= ps_se) ? 1 : 0;
            if (cur_se == best_se) step1.push_back(h);
        }
    }

    // Replicate IRanges::findOverlaps_GNCList + S4Vectors::remapHits tie order:
    //   1) dist ASC
    //   2) first_kq ASC (unlisted_query exon index order)
    //   3) if N_q_exons_chr >= N_s_exons_chr: NCList on subject -> (first_s_start ASC, first_s_end DESC, subject_idx ASC)
    //      else: NCList on query -> subject_idx ASC
    std::sort(step1.begin(), step1.end(), [use_nclist_subj](const DistHit &a, const DistHit &b) {
        if (a.dist != b.dist) return a.dist < b.dist;
        if (a.first_kq != b.first_kq) return a.first_kq < b.first_kq;
        if (use_nclist_subj) {
            if (a.first_s_start != b.first_s_start) return a.first_s_start < b.first_s_start;
            if (a.first_s_end != b.first_s_end) return a.first_s_end > b.first_s_end;
        }
        return a.subject_idx < b.subject_idx;
    });
    hits.swap(step1);
}

} // namespace

void calculate_dist_to_annotation_core(
    const FlatCSRAnnotation &ann,
    const std::vector<std::string> &q_seq_names,
    const std::vector<int32_t> &q_tx_ex_ptr,
    const std::vector<int32_t> &q_ex_seq_0b,
    const std::vector<int32_t> &q_ex_start,
    const std::vector<int32_t> &q_ex_end,
    const std::vector<int8_t> &q_ex_strand,
    const std::vector<int32_t> &q_ex_rank,
    const std::vector<int32_t> &q_ex_end_rank,
    int32_t max_dist,
    int32_t primary_secondary_dist,
    int32_t primary_secondary_dist_se,
    bool ignore_strand,
    int n_threads,
    std::vector<std::vector<DistHit>> &out_per_query
) {
    int32_t n_q = static_cast<int32_t>(q_tx_ex_ptr.size()) - 1;
    out_per_query.assign(n_q, std::vector<DistHit>{});
    if (n_q <= 0 || ann.n_tx == 0) return;

    std::vector<int32_t> q_to_csr_seq(q_seq_names.size(), -1);
    for (size_t i = 0; i < q_seq_names.size(); ++i) {
        auto it = ann.seq_to_id.find(q_seq_names[i]);
        if (it != ann.seq_to_id.end()) q_to_csr_seq[i] = it->second;
    }

    int32_t n_ann_seq = static_cast<int32_t>(ann.seq_names.size());
    std::vector<int64_t> s_ex_cnt(n_ann_seq, 0);
    for (int32_t sid = 0; sid < ann.n_tx; ++sid) {
        int32_t csid = ann.tx_seqid[sid];
        if (csid >= 0 && csid < n_ann_seq) {
            s_ex_cnt[csid] += (ann.tx_exon_ptr[sid + 1] - ann.tx_exon_ptr[sid]);
        }
    }

    struct QueryPrepared {
        int32_t csid = -1;
        int8_t qst = 0;
        int32_t clen = -1;
        int32_t q_st_s = 0, q_st_e = 0, q_en_s = 0, q_en_e = 0;
        int32_t q_cut_min = 0, q_cut_max = 0, q_full_min = 0, q_full_max = 0;
        std::vector<Interval2> q_cut_raw, q_cut_merged, q_ov, q_splice;
    };
    std::vector<QueryPrepared> qprep(n_q);
    std::vector<int64_t> step1_q_cnt(n_ann_seq, 0);

    int act_thr = std::min(std::max(1, n_threads), 16);
    if (n_q < 256) act_thr = 1;

    auto run_parallel = [&](auto &&fn) {
        if (act_thr <= 1) {
            fn(0, n_q, 0);
            return;
        }
        std::vector<std::thread> workers;
        workers.reserve(act_thr);
        for (int t = 0; t < act_thr; ++t) {
            int32_t beg = static_cast<int32_t>((static_cast<int64_t>(n_q) * t) / act_thr);
            int32_t end = static_cast<int32_t>((static_cast<int64_t>(n_q) * (t + 1)) / act_thr);
            workers.emplace_back([&, beg, end, t]() { fn(beg, end, t); });
        }
        for (auto &th : workers) th.join();
    };

    // Precompute query cut/merged/splice geometry in parallel
    std::vector<std::vector<int64_t>> thr_q_cnt(act_thr, std::vector<int64_t>(n_ann_seq, 0));
    run_parallel([&](int32_t beg, int32_t end, int t) {
        auto &loc_cnt = thr_q_cnt[t];
        for (int32_t q = beg; q < end; ++q) {
            int32_t qb = q_tx_ex_ptr[q];
            int32_t qe = q_tx_ex_ptr[q + 1];
            if (qb >= qe) continue;
            int32_t csid = q_to_csr_seq[q_ex_seq_0b[qb]];
            if (csid < 0 || csid >= n_ann_seq) continue;

            auto &qp = qprep[q];
            qp.csid = csid;
            qp.qst = q_ex_strand[qb];
            qp.clen = (csid < static_cast<int32_t>(ann.seq_lengths.size())) ? ann.seq_lengths[csid] : -1;

            int32_t qn = qe - qb;
            qp.q_cut_raw.resize(qn);
            qp.q_st_s = q_ex_start[qb]; qp.q_st_e = q_ex_end[qb];
            qp.q_en_s = q_ex_start[qb]; qp.q_en_e = q_ex_end[qb];
            qp.q_full_min = q_ex_start[qb]; qp.q_full_max = q_ex_end[qb];
            qp.q_ov.clear();

            for (int32_t k = qb; k < qe; ++k) {
                int32_t idx = k - qb;
                int32_t s0 = q_ex_start[k];
                int32_t e0 = q_ex_end[k];
                if (s0 < qp.q_full_min) qp.q_full_min = s0;
                if (e0 > qp.q_full_max) qp.q_full_max = e0;
                if (q_ex_rank[k] == 1) { qp.q_st_s = s0; qp.q_st_e = e0; }
                if (q_ex_end_rank[k] == 1) { qp.q_en_s = s0; qp.q_en_e = e0; }

                int32_t orig_w = e0 - s0 + 1;
                int32_t cs = s0, ce = e0;
                bool is_start = (qp.qst != -1) ? (q_ex_rank[k] == 1) : (q_ex_end_rank[k] == 1);
                if (is_start) cs = ce - 1;
                bool is_end = (qp.qst == -1) ? (q_ex_rank[k] == 1) : (q_ex_end_rank[k] == 1);
                if (is_end) ce = cs + 1;
                qp.q_cut_raw[idx] = Interval2{cs, ce};
                if (orig_w >= max_dist) qp.q_ov.push_back(Interval2{cs, ce});
            }
            loc_cnt[csid] += static_cast<int64_t>(qp.q_ov.size());

            qp.q_cut_merged = qp.q_cut_raw;
            merge_disjoint(qp.q_cut_merged);
            qp.q_splice.clear();
            for (size_t i = 0; i + 1 < qp.q_cut_merged.size(); ++i) {
                int32_t gs = qp.q_cut_merged[i].end + 1;
                int32_t ge = qp.q_cut_merged[i + 1].start - 1;
                if (gs <= ge) qp.q_splice.push_back(Interval2{gs, ge});
            }
            qp.q_cut_min = qp.q_cut_merged.front().start;
            qp.q_cut_max = qp.q_cut_merged.back().end;
        }
    });
    for (int t = 0; t < act_thr; ++t) {
        for (int32_t c = 0; c < n_ann_seq; ++c) step1_q_cnt[c] += thr_q_cnt[t][c];
    }
    std::vector<uint8_t> step1_subj_nclist(n_ann_seq, 0);
    for (int32_t c = 0; c < n_ann_seq; ++c) {
        step1_subj_nclist[c] = (step1_q_cnt[c] >= s_ex_cnt[c]) ? 1 : 0;
    }

    // ---------- Step 1: drop < maxDist, cutStartEnd, within vs subjectExtend(+maxDist) ----------
    run_parallel([&](int32_t beg, int32_t end, int) {
        std::vector<int32_t> cands;
        std::vector<DistHit> hits;
        cands.reserve(64);
        hits.reserve(32);
        for (int32_t q = beg; q < end; ++q) {
            const auto &qp = qprep[q];
            if (qp.csid < 0 || qp.q_ov.empty()) continue;
            hits.clear();
            int32_t ov_min = qp.q_ov[0].start, ov_max = qp.q_ov[0].end;
            for (const auto &r : qp.q_ov) {
                if (r.start < ov_min) ov_min = r.start;
                if (r.end > ov_max) ov_max = r.end;
            }
            cands.clear();
            ann.query_buckets(
                ann.span_index, qp.csid, qp.qst, ignore_strand,
                std::max(1, ov_min - max_dist), ov_max + max_dist,
                [&](int32_t tx) { cands.push_back(tx); });
            std::sort(cands.begin(), cands.end());
            cands.erase(std::unique(cands.begin(), cands.end()), cands.end());

            for (int32_t sid : cands) {
                int32_t sb = ann.tx_exon_ptr[sid], se = ann.tx_exon_ptr[sid + 1];
                bool all_within = true;
                int32_t best_s0 = 0, best_e0 = 0;
                bool has_k0 = false;
                for (size_t kq = 0; kq < qp.q_ov.size(); ++kq) {
                    const auto &qe_iv = qp.q_ov[kq];
                    bool ok_one = false;
                    for (int32_t k = sb; k < se; ++k) {
                        int32_t ex_s = std::max(1, ann.exon_start[k] - max_dist);
                        int32_t ex_e = (qp.clen > 0) ? std::min(qp.clen, ann.exon_end[k] + max_dist)
                                                     : (ann.exon_end[k] + max_dist);
                        if (ex_s <= qe_iv.start && qe_iv.end <= ex_e) {
                            ok_one = true;
                            if (kq == 0) {
                                if (!has_k0 || ex_s < best_s0 || (ex_s == best_s0 && ex_e > best_e0)) {
                                    best_s0 = ex_s;
                                    best_e0 = ex_e;
                                    has_k0 = true;
                                }
                            } else {
                                break;
                            }
                        }
                    }
                    if (!ok_one) { all_within = false; break; }
                }
                if (!all_within) continue;
                DistHit dh;
                score_single_hit(ann, sid, qp.qst, qp.q_cut_merged, qp.q_cut_raw, qp.q_splice,
                                 qp.q_st_s, qp.q_st_e, qp.q_en_s, qp.q_en_e, max_dist, false, dh);
                dh.step = 1;
                dh.first_kq = 0;
                dh.first_s_start = best_s0;
                dh.first_s_end = best_e0;
                hits.push_back(dh);
            }
            if (!hits.empty()) {
                filter_hits_for_query(hits, primary_secondary_dist, primary_secondary_dist_se,
                                      false, step1_subj_nclist[qp.csid] != 0);
                out_per_query[q] = std::move(hits);
            }
        }
    });

    // Count Step 2 remaining queries per chromosome
    std::vector<int64_t> step2_q_cnt(n_ann_seq, 0);
    for (int32_t q = 0; q < n_q; ++q) {
        if (out_per_query[q].empty() && qprep[q].csid >= 0) {
            step2_q_cnt[qprep[q].csid] += (q_tx_ex_ptr[q + 1] - q_tx_ex_ptr[q]);
        }
    }
    std::vector<uint8_t> step2_subj_nclist(n_ann_seq, 0);
    for (int32_t c = 0; c < n_ann_seq; ++c) {
        step2_subj_nclist[c] = (step2_q_cnt[c] >= s_ex_cnt[c]) ? 1 : 0;
    }

    // ---------- Step 2: cutStartEnd, maxDist=0, type="any" ----------
    run_parallel([&](int32_t beg, int32_t end, int) {
        std::vector<int32_t> cands;
        std::vector<DistHit> hits;
        cands.reserve(64);
        hits.reserve(32);
        for (int32_t q = beg; q < end; ++q) {
            if (!out_per_query[q].empty()) continue;
            const auto &qp = qprep[q];
            if (qp.csid < 0) continue;
            hits.clear();
            cands.clear();
            ann.query_buckets(
                ann.span_index, qp.csid, qp.qst, ignore_strand, qp.q_cut_min, qp.q_cut_max,
                [&](int32_t tx) { cands.push_back(tx); });
            std::sort(cands.begin(), cands.end());
            cands.erase(std::unique(cands.begin(), cands.end()), cands.end());
            for (int32_t sid : cands) {
                int32_t sb = ann.tx_exon_ptr[sid], se = ann.tx_exon_ptr[sid + 1];
                int32_t first_kq = -1, best_s = 0, best_e = 0;
                for (size_t kq = 0; kq < qp.q_cut_raw.size(); ++kq) {
                    const auto &qc = qp.q_cut_raw[kq];
                    bool matched_kq = false;
                    for (int32_t k = sb; k < se; ++k) {
                        int32_t ex_s = std::max(1, ann.exon_start[k]);
                        int32_t ex_e = (qp.clen > 0) ? std::min(qp.clen, ann.exon_end[k]) : ann.exon_end[k];
                        if (ex_s <= qc.end && qc.start <= ex_e) {
                            if (!matched_kq || ex_s < best_s || (ex_s == best_s && ex_e > best_e)) {
                                best_s = ex_s;
                                best_e = ex_e;
                            }
                            matched_kq = true;
                        }
                    }
                    if (matched_kq) {
                        first_kq = static_cast<int32_t>(kq);
                        break;
                    }
                }
                if (first_kq < 0) continue;
                DistHit dh;
                score_single_hit(ann, sid, qp.qst, qp.q_cut_merged, qp.q_cut_raw, qp.q_splice,
                                 qp.q_st_s, qp.q_st_e, qp.q_en_s, qp.q_en_e, 0, false, dh);
                dh.step = 2;
                dh.first_kq = first_kq;
                dh.first_s_start = best_s;
                dh.first_s_end = best_e;
                hits.push_back(dh);
            }
            if (!hits.empty()) {
                filter_hits_for_query(hits, primary_secondary_dist, primary_secondary_dist_se,
                                      false, step2_subj_nclist[qp.csid] != 0);
                out_per_query[q] = std::move(hits);
            }
        }
    });

    // Count Step 3 remaining queries per chromosome
    std::vector<int64_t> step3_q_cnt(n_ann_seq, 0);
    for (int32_t q = 0; q < n_q; ++q) {
        if (out_per_query[q].empty() && qprep[q].csid >= 0) {
            step3_q_cnt[qprep[q].csid] += (q_tx_ex_ptr[q + 1] - q_tx_ex_ptr[q]);
        }
    }
    std::vector<uint8_t> step3_subj_nclist(n_ann_seq, 0);
    for (int32_t c = 0; c < n_ann_seq; ++c) {
        step3_subj_nclist[c] = (step3_q_cnt[c] >= s_ex_cnt[c]) ? 1 : 0;
    }

    // ---------- Step 3: no cutStartEnd for overlap, maxDist=0, type="any", DistCalculated=TRUE ----------
    run_parallel([&](int32_t beg, int32_t end, int) {
        std::vector<int32_t> cands;
        std::vector<DistHit> hits;
        cands.reserve(64);
        hits.reserve(32);
        for (int32_t q = beg; q < end; ++q) {
            if (!out_per_query[q].empty()) continue;
            const auto &qp = qprep[q];
            if (qp.csid < 0) continue;
            int32_t qb = q_tx_ex_ptr[q], qe = q_tx_ex_ptr[q + 1];
            hits.clear();
            cands.clear();
            ann.query_buckets(
                ann.span_index, qp.csid, qp.qst, ignore_strand, qp.q_full_min, qp.q_full_max,
                [&](int32_t tx) { cands.push_back(tx); });
            std::sort(cands.begin(), cands.end());
            cands.erase(std::unique(cands.begin(), cands.end()), cands.end());
            for (int32_t sid : cands) {
                int32_t sb = ann.tx_exon_ptr[sid], se = ann.tx_exon_ptr[sid + 1];
                int32_t first_kq = -1, best_s = 0, best_e = 0;
                for (int32_t kq = qb; kq < qe; ++kq) {
                    bool matched_kq = false;
                    for (int32_t k = sb; k < se; ++k) {
                        int32_t ex_s = std::max(1, ann.exon_start[k]);
                        int32_t ex_e = (qp.clen > 0) ? std::min(qp.clen, ann.exon_end[k]) : ann.exon_end[k];
                        if (ex_s <= q_ex_end[kq] && q_ex_start[kq] <= ex_e) {
                            if (!matched_kq || ex_s < best_s || (ex_s == best_s && ex_e > best_e)) {
                                best_s = ex_s;
                                best_e = ex_e;
                            }
                            matched_kq = true;
                        }
                    }
                    if (matched_kq) {
                        first_kq = kq - qb;
                        break;
                    }
                }
                if (first_kq < 0) continue;
                DistHit dh;
                score_single_hit(ann, sid, qp.qst, qp.q_cut_merged, qp.q_cut_raw, qp.q_splice,
                                 qp.q_st_s, qp.q_st_e, qp.q_en_s, qp.q_en_e, 0, true, dh);
                dh.step = 3;
                dh.first_kq = first_kq;
                dh.first_s_start = best_s;
                dh.first_s_end = best_e;
                hits.push_back(dh);
            }
            if (!hits.empty()) {
                filter_hits_for_query(hits, primary_secondary_dist, primary_secondary_dist_se,
                                      true, step3_subj_nclist[qp.csid] != 0);
                out_per_query[q] = std::move(hits);
            }
        }
    });
}

} // namespace bambu_next

// Exported C entrypoint for Stage 2 / calculateDistToAnnotation
extern "C" SEXP bambu_next_calc_dist_to_annotation_cpp(
    SEXP q_seq_names_sexp,
    SEXP q_tx_ex_counts_sexp,
    SEXP q_ex_seq_1b_sexp,
    SEXP q_ex_start_sexp,
    SEXP q_ex_end_sexp,
    SEXP q_ex_strand_sexp,
    SEXP q_ex_rank_sexp,
    SEXP q_ex_end_rank_sexp,
    SEXP csr_ptr_sexp,
    SEXP max_dist_sexp,
    SEXP ps_dist_sexp,
    SEXP ps_dist_se_sexp,
    SEXP ignore_strand_sexp,
    SEXP n_threads_sexp
) {
    BEGIN_RCPP
    if (Rf_isNull(csr_ptr_sexp)) {
        return Rcpp::List::create(
            Rcpp::Named("queryHits") = Rcpp::IntegerVector(0),
            Rcpp::Named("subjectHits") = Rcpp::IntegerVector(0),
            Rcpp::Named("compatible") = Rcpp::LogicalVector(0),
            Rcpp::Named("equal") = Rcpp::LogicalVector(0),
            Rcpp::Named("dist") = Rcpp::IntegerVector(0),
            Rcpp::Named("startMatch") = Rcpp::LogicalVector(0),
            Rcpp::Named("endMatch") = Rcpp::LogicalVector(0),
            Rcpp::Named("uniqueStartLengthQuery") = Rcpp::IntegerVector(0),
            Rcpp::Named("uniqueStartLengthSubject") = Rcpp::IntegerVector(0),
            Rcpp::Named("uniqueEndLengthQuery") = Rcpp::IntegerVector(0),
            Rcpp::Named("uniqueEndLengthSubject") = Rcpp::IntegerVector(0)
        );
    }
    Rcpp::XPtr<bambu_next::FlatCSRAnnotation> ann(csr_ptr_sexp);
    Rcpp::CharacterVector q_snames_r(q_seq_names_sexp);
    Rcpp::IntegerVector q_cnts_r(q_tx_ex_counts_sexp);
    Rcpp::IntegerVector q_es_r(q_ex_seq_1b_sexp);
    Rcpp::IntegerVector q_s_r(q_ex_start_sexp);
    Rcpp::IntegerVector q_e_r(q_ex_end_sexp);
    Rcpp::IntegerVector q_st_r(q_ex_strand_sexp);
    Rcpp::IntegerVector q_rk_r(q_ex_rank_sexp);
    Rcpp::IntegerVector q_erk_r(q_ex_end_rank_sexp);

    int32_t n_q = static_cast<int32_t>(q_cnts_r.size());
    std::vector<std::string> q_snames(q_snames_r.size());
    for (R_xlen_t i = 0; i < q_snames_r.size(); ++i) q_snames[i] = Rcpp::as<std::string>(q_snames_r[i]);

    std::vector<int32_t> q_ptr(n_q + 1);
    int32_t off = 0;
    for (int32_t i = 0; i < n_q; ++i) {
        q_ptr[i] = off;
        off += q_cnts_r[i];
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

    int32_t max_dist = Rf_asInteger(max_dist_sexp);
    int32_t ps_dist = Rf_asInteger(ps_dist_sexp);
    int32_t ps_dist_se = Rf_asInteger(ps_dist_se_sexp);
    bool ignore_strand = Rcpp::as<bool>(ignore_strand_sexp);
    int n_threads = std::max(1, Rf_asInteger(n_threads_sexp));

    std::vector<std::vector<bambu_next::DistHit>> per_q;
    bambu_next::calculate_dist_to_annotation_core(
        *ann, q_snames, q_ptr, q_es, q_s, q_e, q_st, q_rk, q_erk,
        max_dist, ps_dist, ps_dist_se, ignore_strand, n_threads, per_q
    );

    size_t tot = 0;
    for (const auto &v : per_q) tot += v.size();
    Rcpp::IntegerVector o_q(tot), o_s(tot), o_d(tot), o_usq(tot), o_uss(tot), o_ueq(tot), o_ues(tot);
    Rcpp::LogicalVector o_comp(tot), o_eq(tot), o_sm(tot), o_em(tot);
    size_t pos = 0;
    for (uint8_t step = 1; step <= 3; ++step) {
        for (int32_t q = 0; q < n_q; ++q) {
            if (per_q[q].empty() || per_q[q][0].step != step) continue;
            for (const auto &h : per_q[q]) {
                o_q[pos] = q + 1;
                o_s[pos] = h.subject_idx + 1;
                o_comp[pos] = h.compatible;
                o_eq[pos] = h.equal;
                o_d[pos] = h.dist;
                o_sm[pos] = h.start_match;
                o_em[pos] = h.end_match;
                o_usq[pos] = h.u_start_q;
                o_uss[pos] = h.u_start_s;
                o_ueq[pos] = h.u_end_q;
                o_ues[pos] = h.u_end_s;
                ++pos;
            }
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("queryHits") = o_q,
        Rcpp::Named("subjectHits") = o_s,
        Rcpp::Named("compatible") = o_comp,
        Rcpp::Named("equal") = o_eq,
        Rcpp::Named("dist") = o_d,
        Rcpp::Named("startMatch") = o_sm,
        Rcpp::Named("endMatch") = o_em,
        Rcpp::Named("uniqueStartLengthQuery") = o_usq,
        Rcpp::Named("uniqueStartLengthSubject") = o_uss,
        Rcpp::Named("uniqueEndLengthQuery") = o_ueq,
        Rcpp::Named("uniqueEndLengthSubject") = o_ues
    );
    END_RCPP
}
