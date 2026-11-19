#ifndef BAMBU_NEXT_STAGE1_TYPES_HPP
#define BAMBU_NEXT_STAGE1_TYPES_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace bambu_next {

struct UnsplicedRead {
    int32_t seq_id;
    int32_t start;
    int32_t end;
    int32_t read_id;
    int8_t pos_strand; // 1 = '+', 0 = '-'
};

struct UniqueJunctionRec {
    int32_t seq_id;
    int32_t start;
    int32_t end;
    int32_t score = 0;
    int32_t plus_score = 0;
    int32_t minus_score = 0;
    int32_t start_score = 0;
    int32_t end_score = 0;
    char motif[6] = {'N','N','-','N','N','\0'};
    int8_t splice_strand = 0;       // 1 = '+', -1 = '-', 0 = '*'
    int8_t corrected_strand = 0;    // after junctionStrandCorrection
    int8_t annotated_junction = 0;
    int8_t annotated_start = 0;
    int8_t annotated_end = 0;
    int32_t merged_hc_id = -1;      // 0-based index into unique_junctions, -1 if NA (self is separate)
    int32_t merged_hc_id_no_na = -1; // self index if merged_hc_id == -1
};

struct SplicedReadTable {
    // Compact CSR over spliced reads:
    // n_reads spliced reads; read_junc_id stores global 0-based unique junction IDs
    int32_t n_reads = 0;
    std::vector<int32_t> seq_id;
    std::vector<int32_t> start;
    std::vector<int32_t> end;
    std::vector<int32_t> read_id;
    std::vector<int8_t> pos_strand;
    std::vector<int32_t> junc_ptr;     // size n_reads + 1
    std::vector<int32_t> junc_ids;     // global index into unique_junctions
};

struct BamSampleData {
    std::string run_name;
    bool is_indexed = false;
    bool track_reads = false;
    bool all_bam_sq_in_ref = true;

    std::vector<std::string> seq_names;     // seqlevels order used during Stage 1 construction
    std::vector<std::string> ref_seq_names; // genomeSequence seqlevels for final SE
    std::vector<int32_t> genome_seq_lengths;

    std::vector<UnsplicedRead> unspliced;
    SplicedReadTable spliced;
    std::vector<UniqueJunctionRec> unique_junctions;

    // Optional read tracking metadata (populated only if track_reads == true)
    std::vector<std::string> tracked_read_names;
    std::vector<int32_t> tracked_read_ids;

    // Exact Stage 1 warnings
    std::vector<std::string> warnings;
    int64_t n_bad_reads = 0;
    int64_t n_valid_reads = 0;
};

struct ReadClassSet {
    int32_t n_rc = 0;
    std::vector<std::string> rc_names;
    std::vector<int32_t> seq_id;       // 0-based into seq_names
    std::vector<int8_t> strand_code;   // 1 = '+', -1 = '-', 0 = '*'
    std::vector<int32_t> read_count;
    std::vector<double> start_sd;      // NA_REAL when read_count == 1
    std::vector<double> end_sd;
    std::vector<int32_t> read_count_pos_strand;
    std::vector<std::string> intron_starts; // comma-separated or empty for unspliced NA
    std::vector<std::string> intron_ends;
    std::vector<int8_t> confidence_code;    // 0=highConfidenceJunctionReads, 1=lowConfidenceJunctionReads, 2=unsplicedWithin, 3=unsplicedNew

    // Flattened exons CSR (size n_rc + 1)
    std::vector<int32_t> rc_exon_ptr;
    std::vector<int32_t> exon_start;
    std::vector<int32_t> exon_end;
    std::vector<int32_t> exon_rank;
    std::vector<int32_t> exon_end_rank;

    // Flattened readIds CSR (size n_rc + 1)
    std::vector<int32_t> rc_readid_ptr;
    std::vector<int32_t> read_ids;
};

} // namespace bambu_next

#endif // BAMBU_NEXT_STAGE1_TYPES_HPP
