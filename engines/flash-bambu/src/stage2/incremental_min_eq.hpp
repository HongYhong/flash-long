#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "common/flat_csr_index.hpp"

namespace bambu_next {

// Check cutStartEnd compatibility identical to GenomicRanges within + myCompatibleTranscription
inline bool check_cut_compatible_csr(const FlatCSRAnnotation &ann, int32_t q, int32_t s) noexcept {
    int32_t qb = ann.tx_exon_ptr[q], qe = ann.tx_exon_ptr[q + 1];
    int32_t sb = ann.tx_exon_ptr[s], se = ann.tx_exon_ptr[s + 1];
    int32_t qn = qe - qb, sn = se - sb;
    if (qn == 0 || sn == 0) return false;
    if (ann.tx_cut_min_start[q] < ann.tx_cut_min_start[s] ||
        ann.tx_cut_max_end[q] > ann.tx_cut_max_end[s]) {
        return false;
    }

    const int32_t *qs = ann.cut_exon_start.data() + qb;
    const int32_t *qend = ann.cut_exon_end.data() + qb;
    const int32_t *ss = ann.cut_exon_start.data() + sb;
    const int32_t *send = ann.cut_exon_end.data() + sb;

    int32_t j_hint = 0;
    for (int32_t i = 0; i < qn; ++i) {
        int32_t a = qs[i], b = qend[i];
        while (j_hint < sn && send[j_hint] < a) ++j_hint;
        bool covered = false;
        for (int32_t j = j_hint; j < sn && ss[j] <= a; ++j) {
            if (b <= send[j]) {
                covered = true;
                break;
            }
        }
        if (!covered) return false;
    }

    if (qn >= 2) {
        int32_t cur_end = qend[0];
        int32_t sj = 0;
        for (int32_t i = 1; i < qn; ++i) {
            if (qs[i] > cur_end + 1) {
                int32_t g_s = cur_end + 1;
                int32_t g_e = qs[i] - 1;
                while (sj < sn && send[sj] < g_s) ++sj;
                if (sj < sn && ss[sj] <= g_e) return false;
            }
            if (qend[i] > cur_end) cur_end = qend[i];
        }
    }
    return true;
}

// Compute minEqClass on ext_ann = [0..n_novel-1, n_novel..n_tx-1] incrementally
// using ref_csr->eq_tx shifted by +n_novel for the suffix, and computing all
// edges touching [0..n_novel-1] in O(n_novel * avg_overlap).
inline void compute_incremental_min_eq_classes(
    FlatCSRAnnotation &ext_ann,
    const FlatCSRAnnotation &ref_csr,
    int32_t n_novel,
    int n_threads
) {
    int32_t n_tx = ext_ann.n_tx;
    int32_t n_ref = ref_csr.n_tx;
    std::vector<std::vector<int32_t>> eq_lists(n_tx);

    // 1. Copy reference <-> reference edges shifted by +n_novel
    for (int32_t r = 0; r < n_ref; ++r) {
        int32_t b = ref_csr.eq_ptr[r];
        int32_t e = ref_csr.eq_ptr[r + 1];
        auto &dst = eq_lists[n_novel + r];
        dst.resize(e - b);
        for (int32_t k = b; k < e; ++k) {
            dst[k - b] = ref_csr.eq_tx[k] + n_novel;
        }
    }

    // 2. Compute novel -> all edges AND ref -> novel reverse edges across novel queries [0..n_novel-1]
    int use_thr = std::max(1, n_threads);
    std::vector<std::vector<std::pair<int32_t, int32_t>>> rev_edges_thr(use_thr);

#ifdef _OPENMP
#pragma omp parallel num_threads(use_thr) if (n_novel > 64 && use_thr > 1)
#endif
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto &rev_buf = rev_edges_thr[tid];

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 32)
#endif
        for (int32_t q = 0; q < n_novel; ++q) {
            int32_t qb = ext_ann.tx_exon_ptr[q], qe = ext_ann.tx_exon_ptr[q + 1];
            if (qb >= qe) continue;
            int32_t cqs = ext_ann.tx_cut_min_start[q];
            int32_t cqe = ext_ann.tx_cut_max_end[q];
            auto &hits = eq_lists[q];
            ext_ann.query_buckets(
                ext_ann.cut_span_index, ext_ann.tx_seqid[q], ext_ann.tx_strand[q], false,
                cqs, cqe, [&](int32_t s) {
                    // Forward: novel q compatible with subject s (novel or ref)
                    if (check_cut_compatible_csr(ext_ann, q, s)) {
                        hits.push_back(s);
                    }
                    // Reverse: ref s >= n_novel compatible with novel q
                    if (s >= n_novel && check_cut_compatible_csr(ext_ann, s, q)) {
                        rev_buf.emplace_back(s, q);
                    }
                });
            std::sort(hits.begin(), hits.end());
            hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
        }
    }

    // Apply reverse edges ref -> novel (and keep ref lists sorted & unique)
    std::vector<uint8_t> dirty_ref(n_tx, 0);
    for (const auto &buf : rev_edges_thr) {
        for (const auto &ed : buf) {
            eq_lists[ed.first].push_back(ed.second);
            dirty_ref[ed.first] = 1;
        }
    }
    for (int32_t r = n_novel; r < n_tx; ++r) {
        if (dirty_ref[r]) {
            auto &v = eq_lists[r];
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        }
    }

    ext_ann.eq_ptr.resize(n_tx + 1);
    int64_t total_edges = 0;
    for (int32_t q = 0; q < n_tx; ++q) {
        ext_ann.eq_ptr[q] = static_cast<int32_t>(total_edges);
        total_edges += static_cast<int64_t>(eq_lists[q].size());
    }
    ext_ann.eq_ptr[n_tx] = static_cast<int32_t>(total_edges);
    ext_ann.eq_tx.resize(total_edges);
    for (int32_t q = 0; q < n_tx; ++q) {
        std::copy(eq_lists[q].begin(), eq_lists[q].end(), ext_ann.eq_tx.begin() + ext_ann.eq_ptr[q]);
    }
}

} // namespace bambu_next
