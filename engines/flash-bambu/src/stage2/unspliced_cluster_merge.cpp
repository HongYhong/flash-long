#include <Rcpp.h>

#include <algorithm>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace bambu_next {
namespace {

static inline int strand_rank_u(int8_t s) {
    if (s == 1) return 1;
    if (s == -1) return 2;
    return 3;
}

struct UnsplicedItem {
    int32_t seq_idx; // factor index in seq_levels (0-based)
    int8_t strand;   // 1='+', -1='-', 0='*'
    int32_t start;
    int32_t end;
    double read_count;
    double gene_read_prop;
    double tx_score;
    double tx_score_nofit;
};

} // namespace
} // namespace bambu_next

extern "C" SEXP bambu_next_merge_unspliced_cpp(
    SEXP seq_levels_sexp,
    SEXP seq_idx_1b_sexp,
    SEXP strand_code_sexp,
    SEXP start_sexp,
    SEXP end_sexp,
    SEXP read_count_sexp,
    SEXP gene_read_prop_sexp,
    SEXP tx_score_sexp,
    SEXP tx_score_nofit_sexp,
    SEXP stranded_sexp,
    SEXP min_read_count_sexp,
    SEXP min_read_frac_sexp,
    SEXP min_tx_single_exon_sexp
) {
    BEGIN_RCPP
    Rcpp::CharacterVector seq_lvls(seq_levels_sexp);
    Rcpp::IntegerVector seq_idx(seq_idx_1b_sexp);
    Rcpp::IntegerVector st_code(strand_code_sexp);
    Rcpp::IntegerVector starts(start_sexp);
    Rcpp::IntegerVector ends(end_sexp);
    Rcpp::NumericVector rc_vec(read_count_sexp);
    Rcpp::NumericVector grp_vec(gene_read_prop_sexp);
    Rcpp::NumericVector txs_vec(tx_score_sexp);
    Rcpp::NumericVector txsnf_vec(tx_score_nofit_sexp);
    bool stranded = Rcpp::as<bool>(stranded_sexp);
    double min_rc = Rf_asReal(min_read_count_sexp);
    double min_rf = Rf_asReal(min_read_frac_sexp);
    double min_tx_se = Rf_asReal(min_tx_single_exon_sexp);

    int32_t n = static_cast<int32_t>(starts.size());
    if (n == 0) {
        return Rcpp::List::create(
            Rcpp::Named("chr_idx") = Rcpp::IntegerVector(0),
            Rcpp::Named("strand_code") = Rcpp::IntegerVector(0),
            Rcpp::Named("start") = Rcpp::IntegerVector(0),
            Rcpp::Named("end") = Rcpp::IntegerVector(0),
            Rcpp::Named("readCount") = Rcpp::NumericVector(0),
            Rcpp::Named("maxTxScore") = Rcpp::NumericVector(0),
            Rcpp::Named("maxTxScore_noFit") = Rcpp::NumericVector(0),
            Rcpp::Named("NSampleReadCount") = Rcpp::IntegerVector(0),
            Rcpp::Named("NSampleReadProp") = Rcpp::IntegerVector(0),
            Rcpp::Named("NSampleTxScore") = Rcpp::IntegerVector(0)
        );
    }

    std::vector<bambu_next::UnsplicedItem> items;
    items.reserve(n);
    for (int32_t i = 0; i < n; ++i) {
        bambu_next::UnsplicedItem it;
        it.seq_idx = seq_idx[i] - 1;
        it.strand = stranded ? static_cast<int8_t>(st_code[i]) : 0; // ignore.strand = !stranded sets reduced strand to '*'
        it.start = starts[i];
        it.end = ends[i];
        it.read_count = rc_vec[i];
        it.gene_read_prop = grp_vec[i];
        it.tx_score = txs_vec[i];
        it.tx_score_nofit = txsnf_vec[i];
        items.push_back(it);
    }

    // Sort in exact GenomicRanges::reduce order: (seq_idx, strand [+, -, *], start, end)
    std::sort(items.begin(), items.end(), [](const auto &a, const auto &b) {
        if (a.seq_idx != b.seq_idx) return a.seq_idx < b.seq_idx;
        if (a.strand != b.strand) return bambu_next::strand_rank_u(a.strand) < bambu_next::strand_rank_u(b.strand);
        if (a.start != b.start) return a.start < b.start;
        return a.end < b.end;
    });

    struct Cluster {
        int32_t seq_idx;
        int8_t strand;
        int32_t start;
        int32_t end;
        double sum_rc;
        double sum_w_txs;
        double sum_w_txsnf;
        bool any_na_txs;
        bool any_na_txsnf;
        bool any_na_grp;
        int32_t ns_rc;
        int32_t ns_rp;
        int32_t ns_txs;
    };
    std::vector<Cluster> clusters;
    clusters.reserve(n / 2 + 16);

    for (const auto &it : items) {
        bool new_c = clusters.empty() ||
                     clusters.back().seq_idx != it.seq_idx ||
                     clusters.back().strand != it.strand ||
                     it.start > clusters.back().end + 1;
        if (new_c) {
            clusters.push_back(Cluster{
                it.seq_idx, it.strand, it.start, it.end,
                0.0, 0.0, 0.0, false, false, false, 0, 0, 0
            });
        } else if (it.end > clusters.back().end) {
            clusters.back().end = it.end;
        }
        auto &c = clusters.back();
        c.sum_rc += it.read_count;
        if (Rcpp::NumericVector::is_na(it.tx_score)) c.any_na_txs = true;
        else c.sum_w_txs += it.tx_score * it.read_count;
        if (Rcpp::NumericVector::is_na(it.tx_score_nofit)) c.any_na_txsnf = true;
        else c.sum_w_txsnf += it.tx_score_nofit * it.read_count;
        if (it.read_count >= min_rc) c.ns_rc += 1;
        if (Rcpp::NumericVector::is_na(it.gene_read_prop)) c.any_na_grp = true;
        else if (it.gene_read_prop >= min_rf) c.ns_rp += 1;
        if (!Rcpp::NumericVector::is_na(it.tx_score) && it.tx_score > min_tx_se) c.ns_txs += 1;
    }

    int32_t m = static_cast<int32_t>(clusters.size());
    Rcpp::IntegerVector out_chr(m), out_st(m), out_s(m), out_e(m), out_ns_rc(m), out_ns_rp(m), out_ns_txs(m);
    Rcpp::NumericVector out_rc(m), out_txs(m), out_txsnf(m);
    for (int32_t i = 0; i < m; ++i) {
        const auto &c = clusters[i];
        out_chr[i] = c.seq_idx + 1;
        out_st[i] = static_cast<int32_t>(c.strand);
        out_s[i] = c.start;
        out_e[i] = c.end;
        out_rc[i] = c.sum_rc;
        out_txs[i] = (c.any_na_txs || c.sum_rc <= 0.0) ? NA_REAL : (c.sum_w_txs / c.sum_rc);
        out_txsnf[i] = (c.any_na_txsnf || c.sum_rc <= 0.0) ? NA_REAL : (c.sum_w_txsnf / c.sum_rc);
        out_ns_rc[i] = c.ns_rc;
        out_ns_rp[i] = c.any_na_grp ? NA_INTEGER : c.ns_rp;
        out_ns_txs[i] = c.any_na_txs ? NA_INTEGER : c.ns_txs;
    }

    return Rcpp::List::create(
        Rcpp::Named("chr_idx") = out_chr,
        Rcpp::Named("strand_code") = out_st,
        Rcpp::Named("start") = out_s,
        Rcpp::Named("end") = out_e,
        Rcpp::Named("readCount") = out_rc,
        Rcpp::Named("maxTxScore") = out_txs,
        Rcpp::Named("maxTxScore_noFit") = out_txsnf,
        Rcpp::Named("NSampleReadCount") = out_ns_rc,
        Rcpp::Named("NSampleReadProp") = out_ns_rp,
        Rcpp::Named("NSampleTxScore") = out_ns_txs
    );
    END_RCPP
}
