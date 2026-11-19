#include <Rcpp.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace bambu_next {
namespace {

class BufferedFileWriter {
public:
    explicit BufferedFileWriter(const std::string &path) : buf_(1 << 22) { // 4 MiB buffer
        fp_ = std::fopen(path.c_str(), "wb");
        if (!fp_) Rcpp::stop("Failed to open output file for writing: %s", path.c_str());
    }
    ~BufferedFileWriter() {
        if (fp_) {
            flush();
            std::fclose(fp_);
            fp_ = nullptr;
        }
    }
    void write(const char *data, size_t len) {
        if (pos_ + len > buf_.size()) {
            flush();
            if (len >= buf_.size()) {
                if (std::fwrite(data, 1, len, fp_) != len) {
                    Rcpp::stop("Write error to buffered output file");
                }
                return;
            }
        }
        std::memcpy(buf_.data() + pos_, data, len);
        pos_ += len;
    }
    void write_str(const std::string &s) { write(s.data(), s.size()); }
    void write_cstr(const char *s) { write(s, std::strlen(s)); }
    void write_char(char c) { write(&c, 1); }
    void write_int(int32_t v) {
        char tmp[16];
        int n = std::snprintf(tmp, sizeof(tmp), "%d", v);
        if (n > 0) write(tmp, static_cast<size_t>(n));
    }
    void close_ok() {
        flush();
        if (std::fclose(fp_) != 0) {
            fp_ = nullptr;
            Rcpp::stop("Failed to flush/close output file");
        }
        fp_ = nullptr;
    }

private:
    void flush() {
        if (pos_ > 0 && fp_) {
            if (std::fwrite(buf_.data(), 1, pos_, fp_) != pos_) {
                Rcpp::stop("Flush write error to output file");
            }
            pos_ = 0;
        }
    }
    FILE *fp_ = nullptr;
    std::vector<char> buf_;
    size_t pos_ = 0;
};

} // namespace
} // namespace bambu_next

extern "C" SEXP bambu_next_write_gtf_cpp(
    SEXP file_path_sexp,
    SEXP seq_names_sexp,
    SEXP tx_names_sexp,
    SEXP tx_gene_names_sexp,
    SEXP tx_ex_counts_sexp,
    SEXP ex_seq_1b_sexp,
    SEXP ex_start_sexp,
    SEXP ex_end_sexp,
    SEXP ex_strand_code_sexp,
    SEXP ex_rank_sexp
) {
    BEGIN_RCPP
    std::string path = Rcpp::as<std::string>(file_path_sexp);
    Rcpp::CharacterVector seq_nms(seq_names_sexp);
    Rcpp::CharacterVector tx_nms(tx_names_sexp);
    Rcpp::CharacterVector gene_nms(tx_gene_names_sexp);
    Rcpp::IntegerVector tx_cnts(tx_ex_counts_sexp);
    Rcpp::IntegerVector ex_seq(ex_seq_1b_sexp);
    Rcpp::IntegerVector ex_start(ex_start_sexp);
    Rcpp::IntegerVector ex_end(ex_end_sexp);
    Rcpp::IntegerVector ex_st(ex_strand_code_sexp);
    Rcpp::IntegerVector ex_rank(ex_rank_sexp);

    int32_t n_tx = static_cast<int32_t>(tx_cnts.size());
    std::vector<std::string> seqs(seq_nms.size());
    for (R_xlen_t i = 0; i < seq_nms.size(); ++i) seqs[i] = Rcpp::as<std::string>(seq_nms[i]);

    // Build rows matching exact writeToGTF() ordering:
    // arrange(as.character(seqnames), start) across (dfTx + dfExon)
    struct GtfLineRef {
        std::string seq;
        int32_t start;
        int32_t end;
        int8_t strand;
        bool is_tx;
        int32_t tx_idx;
        int32_t rank;
    };
    std::vector<GtfLineRef> lines;
    lines.reserve(static_cast<size_t>(ex_start.size()) + static_cast<size_t>(n_tx));

    int32_t off = 0;
    for (int32_t tx = 0; tx < n_tx; ++tx) {
        int32_t ec = tx_cnts[tx];
        if (ec <= 0) continue;
        int32_t mins = ex_start[off], maxe = ex_end[off];
        int32_t tx_seq = ex_seq[off] - 1;
        int8_t tx_st = static_cast<int8_t>(ex_st[off]);
        for (int32_t k = 0; k < ec; ++k) {
            int32_t p = off + k;
            if (ex_start[p] < mins) mins = ex_start[p];
            if (ex_end[p] > maxe) maxe = ex_end[p];
            lines.push_back(GtfLineRef{
                seqs[ex_seq[p] - 1], ex_start[p], ex_end[p],
                static_cast<int8_t>(ex_st[p]), false, tx, ex_rank[p]
            });
        }
        lines.push_back(GtfLineRef{
            seqs[tx_seq], mins, maxe, tx_st, true, tx, 0
        });
        off += ec;
    }

    // Replicate: rbind(dfTx, dfExon) %>% group_by(group_name) %>% arrange(as.character(seqnames), start)
    // Wait: stable_sort on (seq, start) preserves relative order within ties (dfTx before dfExon since order isn't critical? Wait dfTx was FIRST in rbind(dfTx, dfExon)! Put tx BEFORE ex!).
    std::stable_sort(lines.begin(), lines.end(), [](const GtfLineRef &a, const GtfLineRef &b) {
        if (a.seq != b.seq) return a.seq < b.seq;
        if (a.start != b.start) return a.start < b.start;
        return a.is_tx && !b.is_tx;
    });

    std::vector<std::string> tx_s(n_tx), gn_s(n_tx);
    for (int32_t i = 0; i < n_tx; ++i) {
        tx_s[i] = Rcpp::as<std::string>(tx_nms[i]);
        gn_s[i] = Rcpp::as<std::string>(gene_nms[i]);
    }

    bambu_next::BufferedFileWriter w(path);
    for (const auto &ln : lines) {
        w.write_str(ln.seq);
        w.write_cstr("\tBambu\t");
        w.write_cstr(ln.is_tx ? "transcript\t" : "exon\t");
        w.write_int(ln.start);
        w.write_char('\t');
        w.write_int(ln.end);
        w.write_cstr("\t.\t");
        char sc = (ln.strand == 1) ? '+' : ((ln.strand == -1) ? '-' : '.');
        w.write_char(sc);
        w.write_cstr("\t.\tgene_id \"");
        w.write_str(gn_s[ln.tx_idx]);
        w.write_cstr("\"; transcript_id \"");
        w.write_str(tx_s[ln.tx_idx]);
        if (ln.is_tx) {
            w.write_cstr("\";\n");
        } else {
            w.write_cstr("\"; exon_number \"");
            w.write_int(ln.rank);
            w.write_cstr("\";\n");
        }
    }
    w.close_ok();
    return R_NilValue;
    END_RCPP
}
