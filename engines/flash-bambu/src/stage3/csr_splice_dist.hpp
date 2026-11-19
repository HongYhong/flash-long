#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "common/flat_csr_index.hpp"

namespace bambu_next {

struct Interval2 {
    int32_t start;
    int32_t end;
};

struct DistHit {
    int32_t subject_idx; // 0-based index into ann
    int32_t uq;
    int32_t us;
    int32_t q_out;
    int32_t s_out;
    bool compatible;
    bool equal;
    bool start_match;
    bool end_match;
    int32_t u_start_q;
    int32_t u_start_s;
    int32_t u_end_q;
    int32_t u_end_s;
    int32_t dist;
    uint8_t step; // 1=within(maxDist), 2=any(cut), 3=any(full)
    int32_t first_kq = 0;
    int32_t first_s_start = 0;
    int32_t first_s_end = 0;
};

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
);

} // namespace bambu_next
