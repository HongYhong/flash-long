#pragma once
#include <algorithm>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace bambu_next {

inline int8_t strand_str_to_code(const char *s) noexcept {
    if (!s) return 0;
    if (s[0] == '+') return 1;
    if (s[0] == '-') return -1;
    return 0; // '*' or '.'
}

inline int strand_code_to_bucket(int8_t c) noexcept {
    if (c == 1) return 0;
    if (c == -1) return 1;
    return 2;
}

struct ExonInterval {
    int32_t start;
    int32_t end;
};

// Block-max spatial bucket on a fixed (seqid, strand_bucket) ordered by min_start
struct SpatialBucket {
    static constexpr int32_t BLOCK_SIZE = 64;
    std::vector<int32_t> tx_idx;
    std::vector<int32_t> min_start;
    std::vector<int32_t> max_end;
    std::vector<int32_t> block_max_end;

    void build(const std::vector<int32_t> &indices,
               const std::vector<int32_t> &all_starts,
               const std::vector<int32_t> &all_ends) {
        tx_idx = indices;
        std::sort(tx_idx.begin(), tx_idx.end(), [&](int32_t a, int32_t b) {
            if (all_starts[a] != all_starts[b]) return all_starts[a] < all_starts[b];
            if (all_ends[a] != all_ends[b]) return all_ends[a] < all_ends[b];
            return a < b;
        });
        size_t n = tx_idx.size();
        min_start.resize(n);
        max_end.resize(n);
        size_t nb = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
        block_max_end.assign(nb, INT32_MIN);
        for (size_t i = 0; i < n; ++i) {
            int32_t id = tx_idx[i];
            min_start[i] = all_starts[id];
            max_end[i] = all_ends[id];
            size_t b = i / BLOCK_SIZE;
            if (max_end[i] > block_max_end[b]) block_max_end[b] = max_end[i];
        }
    }

    template <typename Fn>
    void query_overlap(int32_t q_start, int32_t q_end, Fn &&cb) const {
        if (tx_idx.empty()) return;
        auto it = std::upper_bound(min_start.begin(), min_start.end(), q_end);
        int32_t limit = static_cast<int32_t>(it - min_start.begin());
        if (limit <= 0) return;
        int32_t last_block = (limit - 1) / BLOCK_SIZE;
        for (int32_t b = 0; b <= last_block; ++b) {
            if (block_max_end[b] < q_start) continue;
            int32_t beg = b * BLOCK_SIZE;
            int32_t end = std::min(limit, beg + BLOCK_SIZE);
            for (int32_t i = beg; i < end; ++i) {
                if (max_end[i] >= q_start) {
                    cb(tx_idx[i]);
                }
            }
        }
    }
};

struct FlatCSRAnnotation {
    int32_t n_tx = 0;
    int32_t n_seq = 0;
    int32_t n_gene = 0;

    std::vector<std::string> seq_names;
    std::unordered_map<std::string, int32_t> seq_to_id;
    std::vector<int32_t> seq_lengths; // <=0 if unknown

    std::vector<std::string> gene_names;
    std::unordered_map<std::string, int32_t> gene_to_id;

    std::vector<std::string> tx_names;
    std::vector<int32_t> tx_seqid;
    std::vector<int8_t> tx_strand; // +1, -1, 0
    std::vector<int32_t> tx_geneid;

    std::vector<int32_t> tx_min_start;
    std::vector<int32_t> tx_max_end;
    std::vector<int32_t> tx_cut_min_start;
    std::vector<int32_t> tx_cut_max_end;

    // CSR exon arrays (per transcript sorted by start)
    std::vector<int32_t> tx_exon_ptr; // size n_tx + 1
    std::vector<int32_t> exon_start;
    std::vector<int32_t> exon_end;
    std::vector<int32_t> exon_rank;
    std::vector<int32_t> exon_end_rank;

    // CSR cutStartEnd exons (sorted order preserved per tx)
    std::vector<int32_t> cut_exon_start;
    std::vector<int32_t> cut_exon_end;

    // Precomputed minEqClass CSR (0-based subject tx indices, sorted per query tx)
    std::vector<int32_t> eq_ptr; // size n_tx + 1
    std::vector<int32_t> eq_tx;

    // Spatial indices: [seqid * 3 + strand_bucket]
    std::vector<SpatialBucket> span_index;
    std::vector<SpatialBucket> cut_span_index;

    void rebuild_derived_indices() {
        cut_exon_start.resize(exon_start.size());
        cut_exon_end.resize(exon_end.size());
        tx_min_start.resize(n_tx);
        tx_max_end.resize(n_tx);
        tx_cut_min_start.resize(n_tx);
        tx_cut_max_end.resize(n_tx);

        for (int32_t tx = 0; tx < n_tx; ++tx) {
            int32_t b = tx_exon_ptr[tx];
            int32_t e = tx_exon_ptr[tx + 1];
            int32_t mins = INT32_MAX, maxe = INT32_MIN;
            int32_t cmins = INT32_MAX, cmaxe = INT32_MIN;
            int32_t n_ex = e - b;
            for (int32_t k = b; k < e; ++k) {
                if (exon_start[k] < mins) mins = exon_start[k];
                if (exon_end[k] > maxe) maxe = exon_end[k];
                cut_exon_start[k] = exon_start[k];
                cut_exon_end[k] = exon_end[k];
            }
            if (n_ex > 0) {
                int8_t st = tx_strand[tx];
                // Replicate cutStartEndFromGrangesList():
                // startExonsSet: (exon_rank == 1 & strand != "-") | (exon_endRank == 1 & strand == "-")
                // endExonsSet:   (exon_rank == 1 & strand == "-") | (exon_endRank == 1 & strand != "-")
                for (int32_t k = b; k < e; ++k) {
                    bool is_start_ex = (st != -1) ? (exon_rank[k] == 1) : (exon_end_rank[k] == 1);
                    if (is_start_ex) {
                        cut_exon_start[k] = cut_exon_end[k] - 1;
                    }
                }
                for (int32_t k = b; k < e; ++k) {
                    bool is_end_ex = (st == -1) ? (exon_rank[k] == 1) : (exon_end_rank[k] == 1);
                    if (is_end_ex) {
                        cut_exon_end[k] = cut_exon_start[k] + 1;
                    }
                }
                for (int32_t k = b; k < e; ++k) {
                    if (cut_exon_start[k] < cmins) cmins = cut_exon_start[k];
                    if (cut_exon_end[k] > cmaxe) cmaxe = cut_exon_end[k];
                }
            } else {
                mins = maxe = cmins = cmaxe = 0;
            }
            tx_min_start[tx] = mins;
            tx_max_end[tx] = maxe;
            tx_cut_min_start[tx] = cmins;
            tx_cut_max_end[tx] = cmaxe;
        }

        int32_t nb = n_seq * 3;
        std::vector<std::vector<int32_t>> groups(nb);
        for (int32_t tx = 0; tx < n_tx; ++tx) {
            int32_t sid = tx_seqid[tx];
            if (sid >= 0 && sid < n_seq) {
                groups[sid * 3 + strand_code_to_bucket(tx_strand[tx])].push_back(tx);
            }
        }
        span_index.resize(nb);
        cut_span_index.resize(nb);
        for (int32_t b = 0; b < nb; ++b) {
            span_index[b].build(groups[b], tx_min_start, tx_max_end);
            cut_span_index[b].build(groups[b], tx_cut_min_start, tx_cut_max_end);
        }
    }

    template <typename Fn>
    void query_buckets(const std::vector<SpatialBucket> &idx, int32_t seqid, int8_t q_strand,
                       bool ignore_strand, int32_t q_start, int32_t q_end, Fn &&cb) const {
        if (seqid < 0 || seqid >= n_seq) return;
        int32_t base = seqid * 3;
        if (ignore_strand || q_strand == 0) {
            idx[base + 0].query_overlap(q_start, q_end, cb);
            idx[base + 1].query_overlap(q_start, q_end, cb);
            idx[base + 2].query_overlap(q_start, q_end, cb);
        } else if (q_strand == 1) {
            idx[base + 0].query_overlap(q_start, q_end, cb);
            idx[base + 2].query_overlap(q_start, q_end, cb);
        } else {
            idx[base + 1].query_overlap(q_start, q_end, cb);
            idx[base + 2].query_overlap(q_start, q_end, cb);
        }
    }
};

} // namespace bambu_next
