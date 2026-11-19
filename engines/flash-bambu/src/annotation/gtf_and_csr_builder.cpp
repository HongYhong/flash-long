#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include "S4Vectors_interface.h"
#include "IRanges_interface.h"
}

#include "common/flat_csr_index.hpp"

namespace bambu_next {

// Compute cutStartEnd exon compatibility identical to:
//   findOverlaps(cut_q, cut_s, type="within") + myCompatibleTranscription
// Assumptions verified: cut_exons per tx are sorted by start;
// checks exon-level within + setdiff(cut_q, cut_s)==empty + intersect(cut_s, gaps(cut_q))==empty.
inline bool check_cut_compatible(const FlatCSRAnnotation &ann, int32_t q, int32_t s) noexcept {
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

    // Step 1: Every cut exon of q must be contained within >= 1 cut exon of s (GenomicRanges within rule)
    // Because both arrays are ordered by start, advance subject index when send[j] < qs[i]
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

    // Step 2: intersect(cut_s, myGaps(cut_q)) == empty
    // Note myGaps merges overlapping/abutting ranges of q first (like IRanges::gaps)
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

void compute_min_eq_classes(FlatCSRAnnotation &ann, int n_threads) {
    int32_t n_tx = ann.n_tx;
    std::vector<std::vector<int32_t>> eq_lists(n_tx);

    // Parallel loop over transcripts reading immutable span_index/cut arrays
#ifdef _OPENMP
    int use_thr = std::max(1, n_threads);
#pragma omp parallel for schedule(dynamic, 256) num_threads(use_thr) if (n_tx > 1024 && use_thr > 1)
#endif
    for (int32_t q = 0; q < n_tx; ++q) {
        int32_t qb = ann.tx_exon_ptr[q], qe = ann.tx_exon_ptr[q + 1];
        if (qb >= qe) continue;
        int32_t cqs = ann.tx_cut_min_start[q];
        int32_t cqe = ann.tx_cut_max_end[q];
        auto &hits = eq_lists[q];
        ann.query_buckets(
            ann.cut_span_index, ann.tx_seqid[q], ann.tx_strand[q], false,
            cqs, cqe, [&](int32_t s) {
                if (check_cut_compatible(ann, q, s)) {
                    hits.push_back(s);
                }
            });
        std::sort(hits.begin(), hits.end());
        hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
    }

    ann.eq_ptr.resize(n_tx + 1);
    int64_t total_edges = 0;
    for (int32_t q = 0; q < n_tx; ++q) {
        ann.eq_ptr[q] = static_cast<int32_t>(total_edges);
        total_edges += static_cast<int64_t>(eq_lists[q].size());
    }
    ann.eq_ptr[n_tx] = static_cast<int32_t>(total_edges);
    ann.eq_tx.resize(total_edges);
    for (int32_t q = 0; q < n_tx; ++q) {
        std::copy(eq_lists[q].begin(), eq_lists[q].end(), ann.eq_tx.begin() + ann.eq_ptr[q]);
    }
}

namespace {

struct RawExonRec {
    int32_t seq_idx;
    int32_t gene_idx;
    int32_t tx_idx;
    int32_t start;
    int32_t end;
    int8_t strand;
};

// Replicate utils::read.delim(..., quote = "\"") attribute stripping + gsub(".*key (.*?);.*", "\\1", attr)
static inline std::string extract_gtf_attr_bambu(const char *attr_raw, size_t len, const char *key, size_t key_len) {
    // Build stripped string (ignoring '"') OR scan while skipping '"' directly!
    // Scanning directly avoids per-line allocation, matching Bambu's regex .*key ([^;]*);
    std::string cleaned;
    cleaned.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        if (attr_raw[i] != '"') cleaned.push_back(attr_raw[i]);
    }
    // Look for "key " (note Bambu uses ".*gene_id (.*?);.*" matching even substring if greedy .*!)
    // Use greedy last occurrence to match POSIX ERE ".*key (.*?);.*"!
    // Wait: in R gsub('.*gene_id (.*?);.*', '\\1', s): greedy leading .* finds the LAST occurrence if multiple!
    // For standard GTF, gene_id and transcript_id occur at most once per exon row.
    size_t pos = cleaned.rfind(std::string_view(key, key_len));
    while (pos != std::string::npos) {
        if (pos + key_len < cleaned.size() && cleaned[pos + key_len] == ' ') {
            size_t val_start = pos + key_len + 1;
            size_t semi = cleaned.find(';', val_start);
            if (semi != std::string::npos) {
                return cleaned.substr(val_start, semi - val_start);
            }
            return cleaned.substr(val_start);
        }
        if (pos == 0) break;
        pos = cleaned.rfind(std::string_view(key, key_len), pos - 1);
    }
    return cleaned;
}

} // namespace

// Stream parse GTF (supports plain .gtf and gzipped .gtf.gz via htslib htsFile)
static void parse_gtf_to_flat_csr(const std::string &path, FlatCSRAnnotation &ann, int n_threads) {
    htsFile *hf = hts_open(path.c_str(), "r");
    if (!hf) {
        Rcpp::stop("Failed to open GTF file: %s", path.c_str());
    }

    struct TempRow {
        std::string seqname;
        std::string gene_id;
        std::string tx_id;
        int32_t start;
        int32_t end;
        int8_t strand;
    };
    std::vector<TempRow> rows;
    rows.reserve(1 << 16);

    kstring_t ks = {0, 0, nullptr};
    int ret = 0;
    while ((ret = hts_getline(hf, '\n', &ks)) >= 0) {
        if (ks.l == 0 || ks.s[0] == '#') continue;
        // Split up to 9 tab columns
        const char *p = ks.s;
        const char *end = ks.s + ks.l;
        const char *cols[9];
        size_t clen[9];
        int ncol = 0;
        const char *cur = p;
        for (const char *c = p; c <= end && ncol < 9; ++c) {
            if (c == end || *c == '\t') {
                cols[ncol] = cur;
                clen[ncol] = static_cast<size_t>(c - cur);
                ncol++;
                cur = c + 1;
                if (ncol == 8 && c < end) {
                    cols[8] = cur;
                    clen[8] = static_cast<size_t>(end - cur);
                    ncol = 9;
                    break;
                }
            }
        }
        if (ncol < 9) continue;
        // Filter feature == "exon"
        if (clen[2] != 4 || std::memcmp(cols[2], "exon", 4) != 0) continue;

        int32_t s_pos = 0, e_pos = 0;
        {
            char buf[32];
            size_t sl = std::min<size_t>(clen[3], 30);
            std::memcpy(buf, cols[3], sl); buf[sl] = '\0';
            s_pos = std::atoi(buf);
            size_t el = std::min<size_t>(clen[4], 30);
            std::memcpy(buf, cols[4], el); buf[el] = '\0';
            e_pos = std::atoi(buf);
        }
        int8_t st = 0;
        if (clen[6] > 0) {
            if (cols[6][0] == '+') st = 1;
            else if (cols[6][0] == '-') st = -1;
            else st = 0; // '.' -> '*'
        }
        std::string seq(cols[0], clen[0]);
        std::string gid = extract_gtf_attr_bambu(cols[8], clen[8], "gene_id", 7);
        std::string tid = extract_gtf_attr_bambu(cols[8], clen[8], "transcript_id", 13);

        rows.push_back(TempRow{std::move(seq), std::move(gid), std::move(tid), s_pos, e_pos, st});
    }
    if (ks.s) free(ks.s);
    if (hts_close(hf) != 0) {
        Rcpp::warning("Error closing GTF file: %s", path.c_str());
    }

    // Replicate Bambu's multiTxCheck: distinct(seqname, GENEID) -> suffix -1, -2 if same GENEID on >1 seqname
    struct SeqGeneKey {
        std::string seq;
        std::string gene;
        bool operator==(const SeqGeneKey &o) const { return seq == o.seq && gene == o.gene; }
    };
    struct SeqGeneHash {
        size_t operator()(const SeqGeneKey &k) const {
            return std::hash<std::string>{}(k.seq) ^ (std::hash<std::string>{}(k.gene) << 1);
        }
    };
    std::unordered_map<std::string, int32_t> gene_seq_cnt;
    std::unordered_map<SeqGeneKey, int32_t, SeqGeneHash> seq_gene_rank;
    for (const auto &r : rows) {
        SeqGeneKey key{r.seqname, r.gene_id};
        if (seq_gene_rank.find(key) == seq_gene_rank.end()) {
            int32_t idx = ++gene_seq_cnt[r.gene_id];
            seq_gene_rank.emplace(std::move(key), idx);
        }
    }
    bool has_multi_chr_gene = false;
    for (const auto &kv : gene_seq_cnt) {
        if (kv.second > 1) {
            has_multi_chr_gene = true;
            break;
        }
    }
    if (has_multi_chr_gene) {
        Rcpp::warning("Annotations contain duplicated transcript names\n                        Transcript names will be made unique");
        for (auto &r : rows) {
            if (gene_seq_cnt[r.gene_id] > 1) {
                int32_t rk = seq_gene_rank[SeqGeneKey{r.seqname, r.gene_id}];
                std::string suffix = "-" + std::to_string(rk);
                r.gene_id += suffix;
                r.tx_id += suffix;
            }
        }
    }

    // Register unique seq_names in encounter order
    std::unordered_map<std::string, int32_t> tx_map;
    tx_map.reserve(rows.size() / 4 + 16);
    std::vector<std::string> raw_tx_names;
    std::vector<std::string> raw_tx_gene;
    std::vector<int32_t> raw_tx_seq;
    std::vector<int8_t> raw_tx_strand;

    for (const auto &r : rows) {
        int32_t sid;
        auto sit = ann.seq_to_id.find(r.seqname);
        if (sit == ann.seq_to_id.end()) {
            sid = ann.n_seq++;
            ann.seq_to_id.emplace(r.seqname, sid);
            ann.seq_names.push_back(r.seqname);
        } else {
            sid = sit->second;
        }
        auto tit = tx_map.find(r.tx_id);
        if (tit == tx_map.end()) {
            int32_t tidx = static_cast<int32_t>(raw_tx_names.size());
            tx_map.emplace(r.tx_id, tidx);
            raw_tx_names.push_back(r.tx_id);
            raw_tx_gene.push_back(r.gene_id);
            raw_tx_seq.push_back(sid);
            raw_tx_strand.push_back(r.strand);
        }
    }

    int32_t n_tx = static_cast<int32_t>(raw_tx_names.size());
    ann.n_tx = n_tx;

    // Sort transcript indices lexicographically by TXNAME (matching splitAsList character splitting in makeGRangesListFromDataFrame)
    std::vector<int32_t> perm(n_tx);
    for (int32_t i = 0; i < n_tx; ++i) perm[i] = i;
    std::sort(perm.begin(), perm.end(), [&](int32_t a, int32_t b) {
        return raw_tx_names[a] < raw_tx_names[b];
    });
    std::vector<int32_t> old_to_new(n_tx);
    for (int32_t i = 0; i < n_tx; ++i) old_to_new[perm[i]] = i;

    ann.tx_names.resize(n_tx);
    ann.tx_seqid.resize(n_tx);
    ann.tx_strand.resize(n_tx);
    ann.tx_geneid.resize(n_tx);
    std::vector<std::vector<std::pair<int32_t, int32_t>>> tx_exons(n_tx);

    for (int32_t new_id = 0; new_id < n_tx; ++new_id) {
        int32_t old_id = perm[new_id];
        ann.tx_names[new_id] = std::move(raw_tx_names[old_id]);
        ann.tx_seqid[new_id] = raw_tx_seq[old_id];
        ann.tx_strand[new_id] = raw_tx_strand[old_id];
        const std::string &g = raw_tx_gene[old_id];
        auto git = ann.gene_to_id.find(g);
        int32_t gid_idx;
        if (git == ann.gene_to_id.end()) {
            gid_idx = ann.n_gene++;
            ann.gene_to_id.emplace(g, gid_idx);
            ann.gene_names.push_back(g);
        } else {
            gid_idx = git->second;
        }
        ann.tx_geneid[new_id] = gid_idx;
    }

    for (const auto &r : rows) {
        int32_t new_id = old_to_new[tx_map[r.tx_id]];
        tx_exons[new_id].emplace_back(r.start, r.end);
    }

    ann.tx_exon_ptr.resize(n_tx + 1);
    int32_t total_exons = 0;
    for (int32_t tx = 0; tx < n_tx; ++tx) {
        ann.tx_exon_ptr[tx] = total_exons;
        total_exons += static_cast<int32_t>(tx_exons[tx].size());
    }
    ann.tx_exon_ptr[n_tx] = total_exons;

    ann.exon_start.resize(total_exons);
    ann.exon_end.resize(total_exons);
    ann.exon_rank.resize(total_exons);
    ann.exon_end_rank.resize(total_exons);

    for (int32_t tx = 0; tx < n_tx; ++tx) {
        auto &ex = tx_exons[tx];
        std::sort(ex.begin(), ex.end(), [](const auto &a, const auto &b) {
            return a.first < b.first; // stable sort preserves GTF row order on ties like IRanges::order(start(grlist))
        });
        int32_t b = ann.tx_exon_ptr[tx];
        int32_t k = static_cast<int32_t>(ex.size());
        int8_t st = ann.tx_strand[tx];
        for (int32_t i = 0; i < k; ++i) {
            ann.exon_start[b + i] = ex[i].first;
            ann.exon_end[b + i] = ex[i].second;
            if (st == -1) {
                ann.exon_rank[b + i] = k - i;
                ann.exon_end_rank[b + i] = i + 1;
            } else {
                ann.exon_rank[b + i] = i + 1;
                ann.exon_end_rank[b + i] = k - i;
            }
        }
    }

    ann.seq_lengths.assign(ann.n_seq, -1);
    ann.rebuild_derived_indices();
    compute_min_eq_classes(ann, n_threads);
}

} // namespace bambu_next

// R-callable C entrypoints
extern "C" SEXP bambu_next_prepare_gtf_cpp(SEXP gtf_path_sexp, SEXP n_threads_sexp) {
    BEGIN_RCPP
    std::string path = Rcpp::as<std::string>(gtf_path_sexp);
    int n_threads = std::max(1, Rcpp::as<int>(n_threads_sexp));

    auto *ann = new bambu_next::FlatCSRAnnotation();
    bambu_next::parse_gtf_to_flat_csr(path, *ann, n_threads);

    int32_t n_tx = ann->n_tx;
    int32_t total_exons = static_cast<int32_t>(ann->exon_start.size());

    // Build flat vectors to assemble CompressedGRangesList efficiently in R
    Rcpp::IntegerVector ex_start(total_exons);
    Rcpp::IntegerVector ex_end(total_exons);
    Rcpp::IntegerVector ex_rank(total_exons);
    Rcpp::IntegerVector ex_end_rank(total_exons);
    Rcpp::IntegerVector ex_seq_idx(total_exons);
    Rcpp::IntegerVector ex_strand_code(total_exons);
    Rcpp::IntegerVector tx_ex_counts(n_tx);

    for (int32_t tx = 0; tx < n_tx; ++tx) {
        int32_t b = ann->tx_exon_ptr[tx];
        int32_t e = ann->tx_exon_ptr[tx + 1];
        tx_ex_counts[tx] = e - b;
        int32_t seq_1b = ann->tx_seqid[tx] + 1;
        int32_t st = ann->tx_strand[tx];
        for (int32_t k = b; k < e; ++k) {
            ex_start[k] = ann->exon_start[k];
            ex_end[k] = ann->exon_end[k];
            ex_rank[k] = ann->exon_rank[k];
            ex_end_rank[k] = ann->exon_end_rank[k];
            ex_seq_idx[k] = seq_1b;
            ex_strand_code[k] = st;
        }
    }

    Rcpp::CharacterVector seq_names(ann->n_seq);
    for (int32_t i = 0; i < ann->n_seq; ++i) seq_names[i] = ann->seq_names[i];

    Rcpp::CharacterVector tx_names(n_tx);
    Rcpp::CharacterVector tx_gene_names(n_tx);
    for (int32_t i = 0; i < n_tx; ++i) {
        tx_names[i] = ann->tx_names[i];
        tx_gene_names[i] = ann->gene_names[ann->tx_geneid[i]];
    }

    // Export precomputed eqClassById as flat R integer vector + lengths (1-based txid)
    int32_t total_eq = static_cast<int32_t>(ann->eq_tx.size());
    Rcpp::IntegerVector eq_vals(total_eq);
    Rcpp::IntegerVector eq_lens(n_tx);
    for (int32_t i = 0; i < total_eq; ++i) eq_vals[i] = ann->eq_tx[i] + 1;
    for (int32_t tx = 0; tx < n_tx; ++tx) eq_lens[tx] = ann->eq_ptr[tx + 1] - ann->eq_ptr[tx];

    Rcpp::XPtr<bambu_next::FlatCSRAnnotation> ptr(ann, true);

    return Rcpp::List::create(
        Rcpp::Named("seq_names") = seq_names,
        Rcpp::Named("tx_names") = tx_names,
        Rcpp::Named("tx_gene_names") = tx_gene_names,
        Rcpp::Named("tx_ex_counts") = tx_ex_counts,
        Rcpp::Named("ex_start") = ex_start,
        Rcpp::Named("ex_end") = ex_end,
        Rcpp::Named("ex_rank") = ex_rank,
        Rcpp::Named("ex_end_rank") = ex_end_rank,
        Rcpp::Named("ex_seq_idx") = ex_seq_idx,
        Rcpp::Named("ex_strand_code") = ex_strand_code,
        Rcpp::Named("eq_vals") = eq_vals,
        Rcpp::Named("eq_lens") = eq_lens,
        Rcpp::Named("csr_ptr") = ptr
    );
    END_RCPP
}

// Fast builder of FlatCSRAnnotation from an existing prepared GRangesList (e.g. .qs or TxDb result)
// Takes flat pre-extracted integer vectors from R / IRanges C ABI to minimize R dispatch overhead
extern "C" SEXP bambu_next_build_csr_from_flat(
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
    SEXP eq_vals_1b_sexp,
    SEXP eq_lens_sexp,
    SEXP compute_eq_if_null_sexp,
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

    bool compute_eq = Rcpp::as<bool>(compute_eq_if_null_sexp);
    if (!Rf_isNull(eq_vals_1b_sexp) && !Rf_isNull(eq_lens_sexp)) {
        Rcpp::IntegerVector eq_v(eq_vals_1b_sexp);
        Rcpp::IntegerVector eq_l(eq_lens_sexp);
        ann->eq_ptr.resize(n_tx + 1);
        int32_t eq_off = 0;
        for (int32_t i = 0; i < n_tx; ++i) {
            ann->eq_ptr[i] = eq_off;
            eq_off += eq_l[i];
        }
        ann->eq_ptr[n_tx] = eq_off;
        ann->eq_tx.resize(eq_v.size());
        for (R_xlen_t i = 0; i < eq_v.size(); ++i) ann->eq_tx[i] = eq_v[i] - 1;
    } else if (compute_eq) {
        int n_threads = std::max(1, Rcpp::as<int>(n_threads_sexp));
        bambu_next::compute_min_eq_classes(*ann, n_threads);
    }

    Rcpp::XPtr<bambu_next::FlatCSRAnnotation> ptr(ann, true);
    if (Rf_isNull(eq_vals_1b_sexp) && compute_eq) {
        Rcpp::IntegerVector out_vals(ann->eq_tx.size());
        Rcpp::IntegerVector out_lens(n_tx);
        for (size_t i = 0; i < ann->eq_tx.size(); ++i) out_vals[i] = ann->eq_tx[i] + 1;
        for (int32_t i = 0; i < n_tx; ++i) out_lens[i] = ann->eq_ptr[i + 1] - ann->eq_ptr[i];
        return Rcpp::List::create(Rcpp::Named("csr_ptr") = ptr,
                                  Rcpp::Named("eq_vals") = out_vals,
                                  Rcpp::Named("eq_lens") = out_lens);
    }
    return Rcpp::List::create(Rcpp::Named("csr_ptr") = ptr);
    END_RCPP
}

namespace {

struct CsrDigestKey {
    int32_t n_tx = 0;
    int32_t n_ex = 0;
    int32_t eq_tot = -1;
    uint64_t h1 = 0;
    uint64_t h2 = 0;
    std::string first_tx;
    std::string last_tx;

    bool operator==(const CsrDigestKey &o) const {
        return n_tx == o.n_tx && n_ex == o.n_ex && eq_tot == o.eq_tot &&
               h1 == o.h1 && h2 == o.h2 &&
               first_tx == o.first_tx && last_tx == o.last_tx;
    }
};

static inline uint64_t splitmix64_step(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static CsrDigestKey compute_csr_key(
    const int32_t n_tx,
    const int *starts, const int *ends, int32_t n_ex,
    const int *eq_lens, int32_t n_eq_l,
    const int *eq_vals, int32_t n_eq_v,
    const char *first_tx, const char *last_tx) {
    CsrDigestKey k;
    k.n_tx = n_tx;
    k.n_ex = n_ex;
    k.eq_tot = (eq_vals != nullptr) ? n_eq_v : -1;
    k.first_tx = first_tx ? first_tx : "";
    k.last_tx = last_tx ? last_tx : "";

    uint64_t a = splitmix64_step(static_cast<uint64_t>(n_tx) * 1315423911ULL + static_cast<uint64_t>(n_ex));
    uint64_t b = splitmix64_step(static_cast<uint64_t>(n_eq_v) * 2654435761ULL + 0x517cc1b727220a95ULL);
    for (int32_t i = 0; i < n_ex; ++i) {
        uint64_t v = (static_cast<uint64_t>(static_cast<uint32_t>(starts[i])) << 32) |
                     static_cast<uint32_t>(ends[i]);
        a = splitmix64_step(a ^ v);
    }
    if (eq_lens != nullptr) {
        for (int32_t i = 0; i < n_eq_l; ++i) {
            b ^= splitmix64_step(static_cast<uint64_t>(static_cast<uint32_t>(eq_lens[i])) + static_cast<uint64_t>(i));
        }
    }
    if (eq_vals != nullptr) {
        for (int32_t i = 0; i < n_eq_v; ++i) {
            b = splitmix64_step(b ^ (static_cast<uint64_t>(static_cast<uint32_t>(eq_vals[i])) + 0x9e3779b9ULL * i));
        }
    }
    k.h1 = a;
    k.h2 = b;
    return k;
}

struct CsrLruEntry {
    CsrDigestKey key;
    bambu_next::FlatCSRAnnotation *raw_ptr = nullptr;
    SEXP preserved_xptr = R_NilValue;
};

static std::mutex g_csr_cache_mtx;
static std::deque<CsrLruEntry> g_csr_cache;
static constexpr size_t kMaxCsrLruEntries = 4;

static SEXP csr_cache_lookup(const CsrDigestKey &key) {
    std::lock_guard<std::mutex> lk(g_csr_cache_mtx);
    for (auto it = g_csr_cache.rbegin(); it != g_csr_cache.rend(); ++it) {
        if (it->key == key && it->preserved_xptr != R_NilValue &&
            R_ExternalPtrAddr(it->preserved_xptr) != nullptr) {
            return it->preserved_xptr;
        }
    }
    return R_NilValue;
}

static void csr_cache_insert(const CsrDigestKey &key, SEXP xptr) {
    if (TYPEOF(xptr) != EXTPTRSXP || R_ExternalPtrAddr(xptr) == nullptr) return;
    std::lock_guard<std::mutex> lk(g_csr_cache_mtx);
    for (auto it = g_csr_cache.begin(); it != g_csr_cache.end(); ++it) {
        if (it->key == key) {
            if (it->preserved_xptr != xptr) {
                R_ReleaseObject(it->preserved_xptr);
                it->preserved_xptr = xptr;
                R_PreserveObject(xptr);
            }
            return;
        }
    }
    while (g_csr_cache.size() >= kMaxCsrLruEntries) {
        if (g_csr_cache.front().preserved_xptr != R_NilValue) {
            R_ReleaseObject(g_csr_cache.front().preserved_xptr);
        }
        g_csr_cache.pop_front();
    }
    R_PreserveObject(xptr);
    g_csr_cache.push_back(CsrLruEntry{key, nullptr, xptr});
}

} // namespace

extern "C" SEXP bambu_next_csr_is_valid(SEXP ptr_sexp) {
    if (TYPEOF(ptr_sexp) != EXTPTRSXP) return Rf_ScalarLogical(0);
    return Rf_ScalarLogical(R_ExternalPtrAddr(ptr_sexp) != nullptr);
}

extern "C" SEXP bambu_next_csr_cache_lookup(
    SEXP n_tx_sexp, SEXP ex_starts_sexp, SEXP ex_ends_sexp,
    SEXP eq_lens_sexp, SEXP eq_vals_sexp, SEXP first_tx_sexp, SEXP last_tx_sexp) {
    int32_t n_tx = Rf_asInteger(n_tx_sexp);
    int32_t n_ex = static_cast<int32_t>(Rf_length(ex_starts_sexp));
    const int *starts = INTEGER(ex_starts_sexp);
    const int *ends = INTEGER(ex_ends_sexp);
    const int *eq_lens = Rf_isNull(eq_lens_sexp) ? nullptr : INTEGER(eq_lens_sexp);
    int32_t n_eq_l = Rf_isNull(eq_lens_sexp) ? 0 : static_cast<int32_t>(Rf_length(eq_lens_sexp));
    const int *eq_vals = Rf_isNull(eq_vals_sexp) ? nullptr : INTEGER(eq_vals_sexp);
    int32_t n_eq_v = Rf_isNull(eq_vals_sexp) ? 0 : static_cast<int32_t>(Rf_length(eq_vals_sexp));
    const char *ft = Rf_length(first_tx_sexp) > 0 ? CHAR(STRING_ELT(first_tx_sexp, 0)) : "";
    const char *lt = Rf_length(last_tx_sexp) > 0 ? CHAR(STRING_ELT(last_tx_sexp, 0)) : "";
    CsrDigestKey k = compute_csr_key(n_tx, starts, ends, n_ex, eq_lens, n_eq_l, eq_vals, n_eq_v, ft, lt);
    return csr_cache_lookup(k);
}

extern "C" SEXP bambu_next_csr_cache_put(
    SEXP xptr_sexp, SEXP n_tx_sexp, SEXP ex_starts_sexp, SEXP ex_ends_sexp,
    SEXP eq_lens_sexp, SEXP eq_vals_sexp, SEXP first_tx_sexp, SEXP last_tx_sexp) {
    int32_t n_tx = Rf_asInteger(n_tx_sexp);
    int32_t n_ex = static_cast<int32_t>(Rf_length(ex_starts_sexp));
    const int *starts = INTEGER(ex_starts_sexp);
    const int *ends = INTEGER(ex_ends_sexp);
    const int *eq_lens = Rf_isNull(eq_lens_sexp) ? nullptr : INTEGER(eq_lens_sexp);
    int32_t n_eq_l = Rf_isNull(eq_lens_sexp) ? 0 : static_cast<int32_t>(Rf_length(eq_lens_sexp));
    const int *eq_vals = Rf_isNull(eq_vals_sexp) ? nullptr : INTEGER(eq_vals_sexp);
    int32_t n_eq_v = Rf_isNull(eq_vals_sexp) ? 0 : static_cast<int32_t>(Rf_length(eq_vals_sexp));
    const char *ft = Rf_length(first_tx_sexp) > 0 ? CHAR(STRING_ELT(first_tx_sexp, 0)) : "";
    const char *lt = Rf_length(last_tx_sexp) > 0 ? CHAR(STRING_ELT(last_tx_sexp, 0)) : "";
    CsrDigestKey k = compute_csr_key(n_tx, starts, ends, n_ex, eq_lens, n_eq_l, eq_vals, n_eq_v, ft, lt);
    csr_cache_insert(k, xptr_sexp);
    return R_NilValue;
}

extern "C" SEXP bambu_next_csr_cache_clear() {
    std::lock_guard<std::mutex> lk(g_csr_cache_mtx);
    for (auto &e : g_csr_cache) {
        if (e.preserved_xptr != R_NilValue) {
            R_ReleaseObject(e.preserved_xptr);
        }
    }
    g_csr_cache.clear();
    return R_NilValue;
}
