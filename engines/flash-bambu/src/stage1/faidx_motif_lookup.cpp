#include <Rcpp.h>
#include <htslib/faidx.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "stage1/stage1_types.hpp"

namespace bambu_next {

void extract_bam_to_sample_data(
    const std::string &bam_path,
    const std::string &bai_path,
    const std::string &run_name,
    const std::vector<std::string> &genome_seqnames,
    const std::vector<int32_t> &genome_seqlengths,
    const std::vector<std::string> &ann_seqnames,
    bool track_reads,
    int n_threads,
    BamSampleData &out
);

namespace {

// Matches GenomicAlignments::NATURAL_INTRON_MOTIFS = c("GT-AG", "GC-AG", "AT-AC", "AT-AA", "AT-AG")
// and reverse complements c("CT-AC", "CT-GC", "GT-AT", "TT-AT", "CT-AT")
static inline int8_t classify_splice_motif(const char m[6]) {
    if (m[2] != '-') return 0;
    // + strand motifs
    if (std::memcmp(m, "GT-AG", 5) == 0 ||
        std::memcmp(m, "GC-AG", 5) == 0 ||
        std::memcmp(m, "AT-AC", 5) == 0 ||
        std::memcmp(m, "AT-AA", 5) == 0 ||
        std::memcmp(m, "AT-AG", 5) == 0) {
        return 1;
    }
    // - strand motifs (reverse complements)
    if (std::memcmp(m, "CT-AC", 5) == 0 ||
        std::memcmp(m, "CT-GC", 5) == 0 ||
        std::memcmp(m, "GT-AT", 5) == 0 ||
        std::memcmp(m, "TT-AT", 5) == 0 ||
        std::memcmp(m, "CT-AT", 5) == 0) {
        return -1;
    }
    return 0;
}

static inline char up_base(char c) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
}

static inline void fetch_dinuc_fai(faidx_t *fai, const char *chr, int32_t pos1b, char out[2]) {
    out[0] = 'N'; out[1] = 'N';
    if (pos1b < 1) return;
    hts_pos_t len = 0;
    char *s = faidx_fetch_seq64(fai, chr, static_cast<hts_pos_t>(pos1b - 1),
                                static_cast<hts_pos_t>(pos1b), &len);
    if (s != nullptr) {
        if (len >= 1) out[0] = up_base(s[0]);
        if (len >= 2) out[1] = up_base(s[1]);
        free(s);
    }
}

static inline void fetch_dinuc_str(const std::string &seq, int32_t pos1b, char out[2]) {
    out[0] = 'N'; out[1] = 'N';
    if (pos1b < 1 || static_cast<size_t>(pos1b + 1) > seq.size()) return;
    out[0] = up_base(seq[pos1b - 1]);
    out[1] = up_base(seq[pos1b]);
}

} // namespace

void fill_junction_motifs_faidx(
    const std::string &fa_path,
    const std::string &fai_path,
    int n_threads,
    BamSampleData &bsd
) {
    size_t n_j = bsd.unique_junctions.size();
    if (n_j == 0) return;

    int act_thr = std::min(std::max(1, n_threads), 8);
    if (n_j < 4096) act_thr = 1;

    std::vector<std::string> err_msgs(act_thr);
    std::vector<std::thread> workers;
    workers.reserve(act_thr);

    for (int t = 0; t < act_thr; ++t) {
        size_t beg = (n_j * static_cast<size_t>(t)) / static_cast<size_t>(act_thr);
        size_t end = (n_j * static_cast<size_t>(t + 1)) / static_cast<size_t>(act_thr);
        workers.emplace_back([&, t, beg, end]() {
            if (beg >= end) return;
            faidx_t *fai = fai_path.empty()
                ? fai_load3(fa_path.c_str(), nullptr, nullptr, 0)
                : fai_load3(fa_path.c_str(), fai_path.c_str(), nullptr, 0);
            if (!fai) {
                err_msgs[t] = "Failed to load .fai index for FASTA";
                return;
            }
            int32_t prev_seq = -1;
            int32_t prev_start = -1;
            char prev_donor[2] = {'N', 'N'};

            for (size_t i = beg; i < end; ++i) {
                auto &uj = bsd.unique_junctions[i];
                const char *chr = bsd.seq_names[uj.seq_id].c_str();
                if (uj.seq_id != prev_seq || uj.start != prev_start) {
                    fetch_dinuc_fai(fai, chr, uj.start, prev_donor);
                    prev_seq = uj.seq_id;
                    prev_start = uj.start;
                }
                char acc[2];
                fetch_dinuc_fai(fai, chr, uj.end - 1, acc);
                uj.motif[0] = prev_donor[0];
                uj.motif[1] = prev_donor[1];
                uj.motif[2] = '-';
                uj.motif[3] = acc[0];
                uj.motif[4] = acc[1];
                uj.motif[5] = '\0';
                uj.splice_strand = classify_splice_motif(uj.motif);
            }
            fai_destroy(fai);
        });
    }
    for (auto &th : workers) th.join();
    for (const auto &e : err_msgs) {
        if (!e.empty()) Rcpp::stop("FASTA faidx error: %s", e.c_str());
    }
}

void fill_junction_motifs_from_chr_seqs(
    const std::unordered_map<std::string, std::string> &chr_seqs,
    BamSampleData &bsd
) {
    size_t n_j = bsd.unique_junctions.size();
    for (size_t i = 0; i < n_j; ++i) {
        auto &uj = bsd.unique_junctions[i];
        const std::string &chr_name = bsd.seq_names[uj.seq_id];
        auto it = chr_seqs.find(chr_name);
        if (it == chr_seqs.end()) {
            std::memcpy(uj.motif, "NN-NN", 6);
            uj.splice_strand = 0;
            continue;
        }
        char donor[2], acc[2];
        fetch_dinuc_str(it->second, uj.start, donor);
        fetch_dinuc_str(it->second, uj.end - 1, acc);
        uj.motif[0] = donor[0];
        uj.motif[1] = donor[1];
        uj.motif[2] = '-';
        uj.motif[3] = acc[0];
        uj.motif[4] = acc[1];
        uj.motif[5] = '\0';
        uj.splice_strand = classify_splice_motif(uj.motif);
    }
}

} // namespace bambu_next

extern "C" SEXP bambu_next_extract_bam_and_motifs_cpp(
    SEXP bam_path_sexp,
    SEXP bai_path_sexp,
    SEXP run_name_sexp,
    SEXP genome_seqnames_sexp,
    SEXP genome_seqlengths_sexp,
    SEXP ann_seqnames_sexp,
    SEXP fa_path_sexp,
    SEXP fai_path_sexp,
    SEXP track_reads_sexp,
    SEXP n_threads_sexp
) {
    BEGIN_RCPP
    std::string bam_path = Rcpp::as<std::string>(bam_path_sexp);
    std::string bai_path = Rf_isNull(bai_path_sexp) ? "" : Rcpp::as<std::string>(bai_path_sexp);
    std::string run_name = Rcpp::as<std::string>(run_name_sexp);
    std::string fa_path = Rf_isNull(fa_path_sexp) ? "" : Rcpp::as<std::string>(fa_path_sexp);
    std::string fai_path = Rf_isNull(fai_path_sexp) ? "" : Rcpp::as<std::string>(fai_path_sexp);
    bool track_reads = Rcpp::as<bool>(track_reads_sexp);
    int n_threads = std::max(1, Rcpp::as<int>(n_threads_sexp));

    Rcpp::CharacterVector gnames_r(genome_seqnames_sexp);
    Rcpp::IntegerVector glens_r(genome_seqlengths_sexp);
    Rcpp::CharacterVector anames_r(ann_seqnames_sexp);

    std::vector<std::string> gnames(gnames_r.size());
    std::vector<int32_t> glens(glens_r.size());
    for (R_xlen_t i = 0; i < gnames_r.size(); ++i) gnames[i] = Rcpp::as<std::string>(gnames_r[i]);
    for (R_xlen_t i = 0; i < glens_r.size(); ++i) {
        glens[i] = (glens_r[i] == NA_INTEGER) ? -1 : glens_r[i];
    }
    std::vector<std::string> anames(anames_r.size());
    for (R_xlen_t i = 0; i < anames_r.size(); ++i) anames[i] = Rcpp::as<std::string>(anames_r[i]);

    auto *bsd = new bambu_next::BamSampleData();
    bambu_next::extract_bam_to_sample_data(
        bam_path, bai_path, run_name, gnames, glens, anames, track_reads, n_threads, *bsd);

    if (!fa_path.empty()) {
        bambu_next::fill_junction_motifs_faidx(fa_path, fai_path, n_threads, *bsd);
    }

    Rcpp::XPtr<bambu_next::BamSampleData> ptr(bsd, true);
    return Rcpp::List::create(
        Rcpp::Named("bsd_ptr") = ptr,
        Rcpp::Named("n_valid_reads") = static_cast<double>(bsd->n_valid_reads),
        Rcpp::Named("n_spliced") = bsd->spliced.n_reads,
        Rcpp::Named("n_unspliced") = static_cast<int>(bsd->unspliced.size()),
        Rcpp::Named("n_unique_junc") = static_cast<int>(bsd->unique_junctions.size()),
        Rcpp::Named("is_indexed") = bsd->is_indexed,
        Rcpp::Named("warnings") = Rcpp::wrap(bsd->warnings)
    );
    END_RCPP
}

extern "C" SEXP bambu_next_fill_motifs_from_r_cpp(
    SEXP bsd_ptr_sexp,
    SEXP seq_names_sexp,
    SEXP seq_chars_sexp
) {
    BEGIN_RCPP
    Rcpp::XPtr<bambu_next::BamSampleData> bsd(bsd_ptr_sexp);
    Rcpp::CharacterVector nms(seq_names_sexp);
    Rcpp::CharacterVector chrs(seq_chars_sexp);
    std::unordered_map<std::string, std::string> mp;
    mp.reserve(nms.size() * 2 + 8);
    for (R_xlen_t i = 0; i < nms.size(); ++i) {
        mp.emplace(Rcpp::as<std::string>(nms[i]), std::string(CHAR(STRING_ELT(chrs, i))));
    }
    bambu_next::fill_junction_motifs_from_chr_seqs(mp, *bsd);
    return R_NilValue;
    END_RCPP
}

// Export uniqueJunctions table for testing / R compatibility verification
extern "C" SEXP bambu_next_export_unique_junctions_cpp(SEXP bsd_ptr_sexp) {
    BEGIN_RCPP
    Rcpp::XPtr<bambu_next::BamSampleData> bsd(bsd_ptr_sexp);
    int32_t n = static_cast<int32_t>(bsd->unique_junctions.size());
    Rcpp::IntegerVector seq_idx(n), start(n), end(n), score(n), plus_sc(n), minus_sc(n), start_sc(n), end_sc(n);
    Rcpp::IntegerVector sp_st(n), corr_st(n), ann_j(n), ann_s(n), ann_e(n), m_hc(n);
    Rcpp::CharacterVector motif(n);
    for (int32_t i = 0; i < n; ++i) {
        const auto &uj = bsd->unique_junctions[i];
        seq_idx[i] = uj.seq_id + 1;
        start[i] = uj.start;
        end[i] = uj.end;
        score[i] = uj.score;
        plus_sc[i] = uj.plus_score;
        minus_sc[i] = uj.minus_score;
        start_sc[i] = uj.start_score;
        end_sc[i] = uj.end_score;
        sp_st[i] = uj.splice_strand;
        corr_st[i] = uj.corrected_strand;
        ann_j[i] = uj.annotated_junction;
        ann_s[i] = uj.annotated_start;
        ann_e[i] = uj.annotated_end;
        m_hc[i] = (uj.merged_hc_id >= 0) ? (uj.merged_hc_id + 1) : NA_INTEGER;
        motif[i] = uj.motif;
    }
    return Rcpp::List::create(
        Rcpp::Named("seq_names") = Rcpp::wrap(bsd->seq_names),
        Rcpp::Named("seq_idx") = seq_idx,
        Rcpp::Named("start") = start,
        Rcpp::Named("end") = end,
        Rcpp::Named("score") = score,
        Rcpp::Named("plus_score") = plus_sc,
        Rcpp::Named("minus_score") = minus_sc,
        Rcpp::Named("startScore") = start_sc,
        Rcpp::Named("endScore") = end_sc,
        Rcpp::Named("spliceMotif") = motif,
        Rcpp::Named("spliceStrand_code") = sp_st,
        Rcpp::Named("corrected_strand_code") = corr_st,
        Rcpp::Named("annotatedJunction") = ann_j,
        Rcpp::Named("annotatedStart") = ann_s,
        Rcpp::Named("annotatedEnd") = ann_e,
        Rcpp::Named("mergedHighConfJunctionId") = m_hc
    );
    END_RCPP
}
