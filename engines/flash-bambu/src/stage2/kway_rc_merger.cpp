#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace bambu_next {
namespace {

struct SplicedKey {
    std::string chr;
    std::string strand;
    std::string intron_starts;
    std::string intron_ends;
    bool istarts_na = false;
    bool iends_na = false;

    bool operator==(const SplicedKey &o) const {
        return istarts_na == o.istarts_na &&
               iends_na == o.iends_na &&
               chr == o.chr &&
               strand == o.strand &&
               intron_starts == o.intron_starts &&
               intron_ends == o.intron_ends;
    }
};

struct SplicedKeyHash {
    size_t operator()(const SplicedKey &k) const {
        std::hash<std::string> h;
        uint64_t a = h(k.chr) ^ (h(k.strand) * 0x9e3779b97f4a7c15ULL);
        uint64_t b = k.istarts_na ? 0x123456789abcdef0ULL : h(k.intron_starts);
        uint64_t c = k.iends_na ? 0x0fedcba987654321ULL : h(k.intron_ends);
        a ^= (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
        a ^= (c + 0x517cc1b727220a95ULL + (a << 6) + (a >> 2));
        return static_cast<size_t>(a);
    }
};

struct SplicedGroup {
    SplicedKey key;
    std::vector<std::pair<int32_t, int64_t>> start_weights;
    std::vector<std::pair<int32_t, int64_t>> end_weights;
    int64_t total_rc_int = 0;
    double total_rc_dbl = 0.0;
    bool all_rc_int_exact = true;
    int32_t n_sample_rc = 0;
    int32_t n_sample_rp = 0;
    int32_t n_sample_txs = 0;
    bool any_txs_non_na = false;
    double max_tx_score = 0.0;
    bool has_tx_score = false;
    double max_tx_score_nofit = 0.0;
    bool has_tx_score_nofit = false;
};

// Replicate R's:
//   sortVector <- rep(vals, times = as.integer(counts))
//   min(sortVector[sortVector >= quantile(sortVector, probs = 0.5, type = 7)])
// For sorted non-decreasing sortVector of length M >= 1:
//   M = 2r + 1 -> quantile type 7 = x[r] (0-based), first >= x[r] is x[r] == x[M / 2].
//   M = 2r     -> quantile type 7 = 0.5*(x[r-1] + x[r]); if x[r-1]==x[r], first >= q is x[r];
//                 if x[r-1] < x[r], q > x[r-1] so first >= q is also x[r] == x[M / 2]!
static inline int32_t weighted_median_upper_middle(std::vector<std::pair<int32_t, int64_t>> &vw) {
    if (vw.empty()) return NA_INTEGER;
    std::sort(vw.begin(), vw.end(), [](const auto &a, const auto &b) {
        return a.first < b.first;
    });
    int64_t M = 0;
    for (const auto &p : vw) {
        if (p.second > 0) M += p.second;
    }
    if (M <= 0) return NA_INTEGER;
    int64_t target_0b = M / 2;
    int64_t cum = 0;
    for (const auto &p : vw) {
        if (p.second <= 0) continue;
        cum += p.second;
        if (cum > target_0b) return p.first;
    }
    return vw.back().first;
}

static inline std::string sexp_str_or_empty(SEXP s, R_xlen_t idx, bool &is_na) {
    SEXP el = STRING_ELT(s, idx);
    if (el == NA_STRING) {
        is_na = true;
        return std::string();
    }
    is_na = false;
    return std::string(CHAR(el));
}

} // namespace
} // namespace bambu_next

extern "C" SEXP bambu_next_kway_merge_spliced_cpp(
    SEXP samples_list_sexp,
    SEXP min_read_count_sexp,
    SEXP min_read_frac_sexp,
    SEXP min_tx_me_sexp,
    SEXP min_tx_se_sexp
) {
    BEGIN_RCPP
    Rcpp::List samples(samples_list_sexp);
    R_xlen_t n_samples = samples.size();
    double min_rc = Rf_asReal(min_read_count_sexp);
    double min_rf = Rf_asReal(min_read_frac_sexp);
    double min_tx_me = Rf_asReal(min_tx_me_sexp);
    double min_tx_se = Rf_asReal(min_tx_se_sexp);

    std::vector<bambu_next::SplicedGroup> groups;
    std::unordered_map<bambu_next::SplicedKey, int32_t, bambu_next::SplicedKeyHash> key_map;

    for (R_xlen_t s = 0; s < n_samples; ++s) {
        Rcpp::List smp = samples[s];
        Rcpp::IntegerVector ex_cnts = smp["ex_counts"];
        Rcpp::IntegerVector ex_start = smp["ex_start"];
        Rcpp::IntegerVector ex_end = smp["ex_end"];
        Rcpp::CharacterVector chr_vec = smp["chr"];
        Rcpp::CharacterVector strand_vec = smp["strand"];
        Rcpp::CharacterVector is_vec = smp["intronStarts"];
        Rcpp::CharacterVector ie_vec = smp["intronEnds"];
        Rcpp::LogicalVector is_hc = smp["is_hc"];
        Rcpp::NumericVector rc_vec = smp["readCount"];
        Rcpp::NumericVector grp_vec = smp["geneReadProp"];
        Rcpp::NumericVector txs_vec = smp["txScore"];
        Rcpp::NumericVector txsnf_vec = smp["txScore_noFit"];
        Rcpp::IntegerVector numex_vec = smp["numExons"];

        int32_t n_rc = static_cast<int32_t>(ex_cnts.size());
        if (s == 0) {
            groups.reserve(static_cast<size_t>(n_rc) + n_rc / 2 + 64);
            key_map.reserve(static_cast<size_t>(n_rc) * 2 + 64);
        }

        int32_t ex_off = 0;
        for (int32_t i = 0; i < n_rc; ++i) {
            int32_t ec = ex_cnts[i];
            int32_t b = ex_off;
            int32_t e = ex_off + ec;
            ex_off = e;

            if (!is_hc[i] || ec <= 0) continue;
            double rc_val = rc_vec[i];
            if (Rcpp::NumericVector::is_na(rc_val) || rc_val < 1.0) continue;

            int32_t min_s = ex_start[b];
            int32_t max_e = ex_end[b];
            for (int32_t k = b + 1; k < e; ++k) {
                if (ex_start[k] < min_s) min_s = ex_start[k];
                if (ex_end[k] > max_e) max_e = ex_end[k];
            }

            bambu_next::SplicedKey key;
            bool dummy_na = false;
            key.chr = bambu_next::sexp_str_or_empty(chr_vec, i, dummy_na);
            key.strand = bambu_next::sexp_str_or_empty(strand_vec, i, dummy_na);
            key.intron_starts = bambu_next::sexp_str_or_empty(is_vec, i, key.istarts_na);
            key.intron_ends = bambu_next::sexp_str_or_empty(ie_vec, i, key.iends_na);

            int32_t gid = -1;
            auto it = key_map.find(key);
            if (it == key_map.end()) {
                gid = static_cast<int32_t>(groups.size());
                key_map.emplace(key, gid);
                groups.emplace_back();
                groups.back().key = std::move(key);
            } else {
                gid = it->second;
            }
            auto &g = groups[gid];
            int64_t rc_w = static_cast<int64_t>(rc_val);
            g.start_weights.emplace_back(min_s, rc_w);
            g.end_weights.emplace_back(max_e, rc_w);
            g.total_rc_dbl += rc_val;
            if (std::floor(rc_val) == rc_val && rc_val <= 2147483647.0) {
                g.total_rc_int += static_cast<int64_t>(rc_val);
            } else {
                g.all_rc_int_exact = false;
            }
            if (rc_val >= min_rc) g.n_sample_rc += 1;

            double grp_val = grp_vec[i];
            if (!Rcpp::NumericVector::is_na(grp_val) && grp_val >= min_rf) {
                g.n_sample_rp += 1;
            }

            double txs_val = txs_vec[i];
            if (!Rcpp::NumericVector::is_na(txs_val)) {
                g.any_txs_non_na = true;
                if (!g.has_tx_score || txs_val > g.max_tx_score) {
                    g.max_tx_score = txs_val;
                    g.has_tx_score = true;
                }
                int32_t nex = numex_vec[i];
                if ((nex >= 2 && txs_val > min_tx_me) || (nex == 1 && txs_val > min_tx_se)) {
                    g.n_sample_txs += 1;
                }
            }

            double txsnf_val = txsnf_vec[i];
            if (!Rcpp::NumericVector::is_na(txsnf_val)) {
                if (!g.has_tx_score_nofit || txsnf_val > g.max_tx_score_nofit) {
                    g.max_tx_score_nofit = txsnf_val;
                    g.has_tx_score_nofit = true;
                }
            }
        }
    }

    int32_t n_out = static_cast<int32_t>(groups.size());
    Rcpp::IntegerVector out_start(n_out), out_end(n_out);
    Rcpp::NumericVector out_rc(n_out);
    Rcpp::CharacterVector out_is(n_out), out_ie(n_out), out_chr(n_out), out_strand(n_out);
    Rcpp::NumericVector out_max_txs(n_out), out_max_txsnf(n_out);
    Rcpp::IntegerVector out_ns_rc(n_out), out_ns_rp(n_out), out_ns_txs(n_out);

    for (int32_t i = 0; i < n_out; ++i) {
        auto &g = groups[i];
        out_start[i] = bambu_next::weighted_median_upper_middle(g.start_weights);
        out_end[i] = bambu_next::weighted_median_upper_middle(g.end_weights);
        out_rc[i] = g.total_rc_dbl;
        if (g.key.istarts_na) SET_STRING_ELT(out_is, i, NA_STRING);
        else out_is[i] = g.key.intron_starts;
        if (g.key.iends_na) SET_STRING_ELT(out_ie, i, NA_STRING);
        else out_ie[i] = g.key.intron_ends;
        out_chr[i] = g.key.chr;
        out_strand[i] = g.key.strand;
        out_max_txs[i] = g.has_tx_score ? g.max_tx_score : NA_REAL;
        out_max_txsnf[i] = g.has_tx_score_nofit ? g.max_tx_score_nofit : NA_REAL;
        out_ns_rc[i] = g.n_sample_rc;
        out_ns_rp[i] = g.n_sample_rp;
        out_ns_txs[i] = (n_samples == 1 && !g.any_txs_non_na) ? NA_INTEGER : g.n_sample_txs;
    }

    return Rcpp::List::create(
        Rcpp::Named("start") = out_start,
        Rcpp::Named("end") = out_end,
        Rcpp::Named("readCount") = out_rc,
        Rcpp::Named("intronStarts") = out_is,
        Rcpp::Named("intronEnds") = out_ie,
        Rcpp::Named("chr") = out_chr,
        Rcpp::Named("strand") = out_strand,
        Rcpp::Named("maxTxScore") = out_max_txs,
        Rcpp::Named("maxTxScore_noFit") = out_max_txsnf,
        Rcpp::Named("NSampleReadCount") = out_ns_rc,
        Rcpp::Named("NSampleReadProp") = out_ns_rp,
        Rcpp::Named("NSampleTxScore") = out_ns_txs
    );
    END_RCPP
}
