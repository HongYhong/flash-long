#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common/flat_csr_index.hpp"
#include "stage1/stage1_types.hpp"

namespace bambu_next {

void apply_junc_ec_predictions(
    BamSampleData &bsd,
    bool can_run_hc,
    const std::vector<double> &plus_ss, const std::vector<double> &plus_se,
    const std::vector<double> &plus_es, const std::vector<double> &plus_ee,
    const std::vector<double> &minus_ss, const std::vector<double> &minus_se,
    const std::vector<double> &minus_es, const std::vector<double> &minus_ee
);

namespace {

static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static inline int strand_factor_rank(int8_t s) {
    // Factor levels c("+", "-", "*") -> 1, 2, 3
    if (s == 1) return 1;
    if (s == -1) return 2;
    return 3;
}

struct SplicedRcKey {
    int32_t seq_id;
    int8_t strand;
    int8_t conf_code;
    std::vector<int32_t> intron_starts;
    std::vector<int32_t> intron_ends;
    uint64_t prehash;

    bool operator==(const SplicedRcKey &o) const {
        return prehash == o.prehash &&
               seq_id == o.seq_id &&
               strand == o.strand &&
               conf_code == o.conf_code &&
               intron_starts == o.intron_starts &&
               intron_ends == o.intron_ends;
    }
};

struct SplicedRcKeyHash {
    size_t operator()(const SplicedRcKey &k) const {
        return static_cast<size_t>(k.prehash);
    }
};

struct SplicedRcGroup {
    SplicedRcKey key;
    std::vector<int32_t> read_starts;
    std::vector<int32_t> read_ends;
    std::vector<int32_t> read_ids;
    int32_t pos_strand_cnt = 0;

    // Computed summary statistics
    int32_t read_count = 0;
    int32_t rc_start = 0;
    int32_t rc_end = 0;
    double start_sd = NA_REAL;
    double end_sd = NA_REAL;
    std::string intron_starts_str;
    std::string intron_ends_str;
};

static inline double compute_sample_sd(const std::vector<int32_t> &v) {
    size_t n = v.size();
    if (n <= 1) return NA_REAL;
    double sum = 0.0;
    for (int32_t x : v) sum += static_cast<double>(x);
    double mean = sum / static_cast<double>(n);
    double ss = 0.0;
    for (int32_t x : v) {
        double d = static_cast<double>(x) - mean;
        ss += d * d;
    }
    return std::sqrt(ss / static_cast<double>(n - 1));
}

static inline double compute_weighted_sample_sd(
    const std::vector<std::pair<int32_t, int32_t>> &val_cnt,
    int64_t total_n
) {
    if (total_n <= 1) return NA_REAL;
    double sum = 0.0;
    for (const auto &vc : val_cnt) {
        sum += static_cast<double>(vc.first) * static_cast<double>(vc.second);
    }
    double mean = sum / static_cast<double>(total_n);
    double ss = 0.0;
    for (const auto &vc : val_cnt) {
        double d = static_cast<double>(vc.first) - mean;
        ss += (d * d) * static_cast<double>(vc.second);
    }
    return std::sqrt(ss / static_cast<double>(total_n - 1));
}

static inline std::string join_ints_comma(const std::vector<int32_t> &v) {
    std::string s;
    s.reserve(v.size() * 10);
    for (size_t i = 0; i < v.size(); ++i) {
        if (i > 0) s.push_back(',');
        s += std::to_string(v[i]);
    }
    return s;
}

struct RefExon {
    int32_t seq_id;
    int32_t start;
    int32_t end;
    int8_t strand;
};

struct UniqueUnspliced {
    int32_t seq_id;
    int32_t start;
    int32_t end;
    int8_t strand; // 1 for '+', -1 for '-'
    int32_t count;
    std::vector<int32_t> read_ids;
};

struct UnsplicedGroupRec {
    int32_t seq_id;
    int32_t start;
    int32_t end;
    int8_t strand;
    int32_t gid;
    std::string rc_name;
    int32_t read_count = 0;
    int32_t pos_strand_cnt = 0;
    double start_sd = NA_REAL;
    double end_sd = NA_REAL;
    std::vector<std::pair<int32_t, int32_t>> rstart_cnts;
    std::vector<std::pair<int32_t, int32_t>> rend_cnts;
    std::vector<int32_t> read_ids;
};

static void build_unspliced_classes_by_ref(
    const std::vector<UniqueUnspliced> &u_reads,
    const std::vector<int32_t> &active_u_indices,
    const std::vector<RefExon> &refs, // sorted by (seq_id, start, end, strand)
    bool stranded,
    const std::string &conf_name, // "unsplicedWithin" or "unsplicedNew"
    std::vector<uint8_t> &matched_mask,
    std::vector<UnsplicedGroupRec> &out_groups
) {
    out_groups.clear();
    if (active_u_indices.empty() || refs.empty()) return;

    // Sort active queries by (seq_id, start)
    std::vector<int32_t> q_order = active_u_indices;
    std::sort(q_order.begin(), q_order.end(), [&](int32_t a, int32_t b) {
        const auto &ua = u_reads[a];
        const auto &ub = u_reads[b];
        if (ua.seq_id != ub.seq_id) return ua.seq_id < ub.seq_id;
        return ua.start < ub.start;
    });

    // Plane sweep per seq_id
    struct HitAgg {
        int32_t seq_id;
        int8_t eff_strand;
        int32_t max_start;
        int32_t min_end;
        int32_t u_idx;
    };
    std::vector<HitAgg> hit_aggs;
    hit_aggs.reserve(active_u_indices.size());

    size_t qi = 0;
    size_t ri = 0;
    size_t n_q = q_order.size();
    size_t n_r = refs.size();

    while (qi < n_q && ri < n_r) {
        int32_t seq_q = u_reads[q_order[qi]].seq_id;
        int32_t seq_r = refs[ri].seq_id;
        if (seq_q < seq_r) {
            while (qi < n_q && u_reads[q_order[qi]].seq_id == seq_q) qi++;
            continue;
        }
        if (seq_r < seq_q) {
            while (ri < n_r && refs[ri].seq_id == seq_r) ri++;
            continue;
        }
        int32_t seq = seq_q;
        size_t q_beg = qi, r_beg = ri;
        while (qi < n_q && u_reads[q_order[qi]].seq_id == seq) qi++;
        while (ri < n_r && refs[ri].seq_id == seq) ri++;

        // Active set of refs overlapping current start
        // Since refs are sorted by start, advance r_ptr while refs[r_ptr].start <= q.start
        // Active refs pruned when ref.end < q.start
        std::vector<size_t> active_refs;
        active_refs.reserve(256);
        size_t r_ptr = r_beg;

        for (size_t k = q_beg; k < qi; ++k) {
            int32_t u_idx = q_order[k];
            const auto &q = u_reads[u_idx];
            while (r_ptr < ri && refs[r_ptr].start <= q.start) {
                if (refs[r_ptr].end >= q.start) active_refs.push_back(r_ptr);
                r_ptr++;
            }
            // Prune active_refs whose end < q.start
            size_t write_p = 0;
            for (size_t idx = 0; idx < active_refs.size(); ++idx) {
                size_t r_idx = active_refs[idx];
                if (refs[r_idx].end >= q.start) {
                    active_refs[write_p++] = r_idx;
                }
            }
            active_refs.resize(write_p);

            // Check within: ref.start <= q.start && ref.end >= q.end
            // If !stranded, only one eff_strand = 0; if stranded, can match ref.strand == q.strand OR ref.strand == 0
            int32_t m_s[3] = {0, 0, 0};
            int32_t m_e[3] = {0, 0, 0};
            uint8_t has_st[3] = {0, 0, 0};
            for (size_t r_idx : active_refs) {
                const auto &rf = refs[r_idx];
                if (rf.end < q.end) continue;
                int8_t eff = stranded ? rf.strand : 0;
                if (stranded && eff != 0 && eff != q.strand) continue;
                int slot = (eff == 1) ? 0 : ((eff == -1) ? 1 : 2);
                if (!has_st[slot]) {
                    has_st[slot] = 1;
                    m_s[slot] = rf.start;
                    m_e[slot] = rf.end;
                } else {
                    if (rf.start > m_s[slot]) m_s[slot] = rf.start;
                    if (rf.end < m_e[slot]) m_e[slot] = rf.end;
                }
            }
            int8_t slot_strands[3] = {1, -1, 0};
            for (int sl = 0; sl < 3; ++sl) {
                if (has_st[sl]) {
                    matched_mask[u_idx] = 1;
                    hit_aggs.push_back(HitAgg{seq, slot_strands[sl], m_s[sl], m_e[sl], u_idx});
                }
            }
        }
    }

    if (hit_aggs.empty()) return;

    // Group hit_aggs by (seq_id, start, end, eff_strand) matching dplyr::group_by(chr, start, end, strand)
    // Wait: order u_idx in original encounter order (ascending u_idx = original reads.singleExon order)
    // so readIds within each RC are appended in exact R order!
    std::sort(hit_aggs.begin(), hit_aggs.end(), [](const HitAgg &a, const HitAgg &b) {
        if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
        if (a.max_start != b.max_start) return a.max_start < b.max_start;
        if (a.min_end != b.min_end) return a.min_end < b.min_end;
        if (a.eff_strand != b.eff_strand) return strand_factor_rank(a.eff_strand) < strand_factor_rank(b.eff_strand);
        return a.u_idx < b.u_idx;
    });

    int32_t gid = 1;
    size_t i = 0;
    while (i < hit_aggs.size()) {
        size_t s = i;
        const auto &h0 = hit_aggs[i];
        while (i < hit_aggs.size() &&
               hit_aggs[i].seq_id == h0.seq_id &&
               hit_aggs[i].max_start == h0.max_start &&
               hit_aggs[i].min_end == h0.min_end &&
               hit_aggs[i].eff_strand == h0.eff_strand) {
            i++;
        }
        UnsplicedGroupRec g;
        g.seq_id = h0.seq_id;
        g.start = h0.max_start;
        g.end = h0.min_end;
        g.strand = h0.eff_strand;
        g.gid = gid++;
        g.rc_name = "rc" + conf_name + "." + std::to_string(g.gid);

        int64_t tot_cnt = 0;
        int64_t pos_cnt = 0;
        for (size_t k = s; k < i; ++k) {
            const auto &u = u_reads[hit_aggs[k].u_idx];
            tot_cnt += u.count;
            if (u.strand == 1) pos_cnt += u.count;
            g.rstart_cnts.emplace_back(u.start, u.count);
            g.rend_cnts.emplace_back(u.end, u.count);
            g.read_ids.insert(g.read_ids.end(), u.read_ids.begin(), u.read_ids.end());
        }
        g.read_count = static_cast<int32_t>(tot_cnt);
        g.pos_strand_cnt = static_cast<int32_t>(pos_cnt);
        g.start_sd = compute_weighted_sample_sd(g.rstart_cnts, tot_cnt);
        g.end_sd = compute_weighted_sample_sd(g.rend_cnts, tot_cnt);
        std::vector<std::pair<int32_t, int32_t>>().swap(g.rstart_cnts);
        std::vector<std::pair<int32_t, int32_t>>().swap(g.rend_cnts);
        out_groups.push_back(std::move(g));
    }

    // Replicate dplyr::group_by(readClassId) lexicographical sort of rc_name strings!
    std::sort(out_groups.begin(), out_groups.end(), [](const UnsplicedGroupRec &a, const UnsplicedGroupRec &b) {
        return a.rc_name < b.rc_name;
    });
}

} // namespace

// Call 3: Apply junction predictions, run spliced + unspliced ReadClass builder, return full ReadClassSet SE data to R
extern "C" SEXP bambu_next_build_read_classes_cpp(
    SEXP bsd_ptr_sexp,
    SEXP csr_ptr_sexp,
    SEXP stranded_sexp,
    SEXP can_run_hc_sexp,
    SEXP plus_ss_sexp, SEXP plus_se_sexp, SEXP plus_es_sexp, SEXP plus_ee_sexp,
    SEXP minus_ss_sexp, SEXP minus_se_sexp, SEXP minus_es_sexp, SEXP minus_ee_sexp
) {
    BEGIN_RCPP
    Rcpp::XPtr<bambu_next::BamSampleData> bsd(bsd_ptr_sexp);
    const FlatCSRAnnotation *ann = Rf_isNull(csr_ptr_sexp)
        ? nullptr
        : Rcpp::XPtr<FlatCSRAnnotation>(csr_ptr_sexp).get();
    bool stranded = Rcpp::as<bool>(stranded_sexp);
    bool can_run_hc = Rcpp::as<bool>(can_run_hc_sexp);

    if (!bsd->unique_junctions.empty()) {
        auto to_vec = [](SEXP s) -> std::vector<double> {
            if (Rf_isNull(s)) return {};
            Rcpp::NumericVector v(s);
            return std::vector<double>(v.begin(), v.end());
        };
        apply_junc_ec_predictions(
            *bsd, can_run_hc,
            to_vec(plus_ss_sexp), to_vec(plus_se_sexp), to_vec(plus_es_sexp), to_vec(plus_ee_sexp),
            to_vec(minus_ss_sexp), to_vec(minus_se_sexp), to_vec(minus_es_sexp), to_vec(minus_ee_sexp)
        );
    }

    // -------------------------------------------------------------------------
    // 1. Construct Spliced ReadClasses
    // -------------------------------------------------------------------------
    std::vector<SplicedRcGroup> sp_groups;
    if (!bsd->unique_junctions.empty() && bsd->spliced.n_reads > 0) {
        std::unordered_map<SplicedRcKey, int32_t, SplicedRcKeyHash> rc_map;
        rc_map.reserve(bsd->spliced.n_reads / 4 + 256);

        std::vector<int32_t> orig_s, s_buf, e_buf, st_buf;
        std::vector<uint8_t> low_buf, rem_flag;
        orig_s.reserve(64); s_buf.reserve(64); e_buf.reserve(64); st_buf.reserve(64);
        low_buf.reserve(64); rem_flag.reserve(64);

        const auto &uj = bsd->unique_junctions;
        int32_t n_sp = bsd->spliced.n_reads;

        for (int32_t r = 0; r < n_sp; ++r) {
            int32_t jb = bsd->spliced.junc_ptr[r];
            int32_t je = bsd->spliced.junc_ptr[r + 1];
            int32_t m = je - jb;
            if (m <= 0) continue;

            s_buf.resize(m);
            e_buf.resize(m);
            st_buf.resize(m);
            low_buf.resize(m);
            for (int32_t k = 0; k < m; ++k) {
                int32_t j = bsd->spliced.junc_ids[jb + k];
                int32_t c = uj[j].merged_hc_id_no_na;
                s_buf[k] = uj[c].start;
                e_buf[k] = uj[c].end;
                st_buf[k] = uj[c].corrected_strand;
                low_buf[k] = (uj[j].merged_hc_id < 0) ? 1 : 0;
            }

            // Replicate correctIntronRanges microexon removal verbatim
            orig_s = s_buf;
            rem_flag.assign(m, 0);
            for (int32_t k = 0; k + 1 < m; ++k) {
                if (orig_s[k + 1] <= e_buf[k]) {
                    s_buf[k + 1] = orig_s[k];
                    rem_flag[k] = 1;
                }
            }

            int32_t write_k = 0;
            int32_t min_is = INT32_MAX;
            int32_t max_ie = INT32_MIN;
            int32_t plus_j = 0, minus_j = 0;
            bool is_low = false;
            for (int32_t k = 0; k < m; ++k) {
                if (rem_flag[k]) continue;
                s_buf[write_k] = s_buf[k];
                e_buf[write_k] = e_buf[k];
                if (s_buf[k] < min_is) min_is = s_buf[k];
                if (e_buf[k] > max_ie) max_ie = e_buf[k];
                if (st_buf[k] == 1) plus_j++;
                else if (st_buf[k] == -1) minus_j++;
                if (low_buf[k]) is_low = true;
                write_k++;
            }
            s_buf.resize(write_k);
            e_buf.resize(write_k);
            if (write_k == 0) continue;

            int8_t read_st = 0;
            if (!stranded) {
                int32_t diff = minus_j - plus_j;
                if (diff < 0) read_st = 1;
                else if (diff > 0) read_st = -1;
            } else {
                read_st = bsd->spliced.pos_strand[r] ? 1 : -1;
            }
            int8_t conf_c = is_low ? 1 : 0;

            int32_t adj_s = std::min(bsd->spliced.start[r], min_is - 2);
            int32_t adj_e = std::max(bsd->spliced.end[r], max_ie + 2);

            uint64_t h = splitmix64(
                (static_cast<uint64_t>(bsd->spliced.seq_id[r]) << 32) ^
                (static_cast<uint64_t>(static_cast<uint8_t>(read_st + 2)) << 8) ^
                static_cast<uint64_t>(conf_c)
            );
            for (int32_t k = 0; k < write_k; ++k) {
                h = splitmix64(h ^ ((static_cast<uint64_t>(static_cast<uint32_t>(s_buf[k])) << 32) |
                                    static_cast<uint32_t>(e_buf[k])));
            }

            SplicedRcKey key;
            key.seq_id = bsd->spliced.seq_id[r];
            key.strand = read_st;
            key.conf_code = conf_c;
            key.intron_starts = s_buf;
            key.intron_ends = e_buf;
            key.prehash = h;

            auto it = rc_map.find(key);
            int32_t gidx;
            if (it == rc_map.end()) {
                gidx = static_cast<int32_t>(sp_groups.size());
                SplicedRcGroup grp;
                grp.key = std::move(key);
                sp_groups.push_back(std::move(grp));
                rc_map.emplace(sp_groups[gidx].key, gidx);
            } else {
                gidx = it->second;
            }
            auto &g = sp_groups[gidx];
            g.read_starts.push_back(adj_s);
            g.read_ends.push_back(adj_e);
            g.read_ids.push_back(bsd->spliced.read_id[r]);
            if (bsd->spliced.pos_strand[r]) g.pos_strand_cnt++;
        }

        // Compute exact quantile start / end and SD per group
        for (auto &g : sp_groups) {
            int32_t n = static_cast<int32_t>(g.read_starts.size());
            g.read_count = n;
            g.start_sd = compute_sample_sd(g.read_starts);
            g.end_sd = compute_sample_sd(g.read_ends);
            std::sort(g.read_starts.begin(), g.read_starts.end());
            std::sort(g.read_ends.begin(), g.read_ends.end());
            int32_t i20 = std::max(0, std::min(n - 1, static_cast<int32_t>(std::ceil(static_cast<double>(n) / 5.0)) - 1));
            int32_t i80 = std::max(0, std::min(n - 1, static_cast<int32_t>(std::ceil(static_cast<double>(n) / 1.25)) - 1));
            g.rc_start = g.read_starts[i20];
            g.rc_end = g.read_ends[i80];
            g.intron_starts_str = join_ints_comma(g.key.intron_starts);
            g.intron_ends_str = join_ints_comma(g.key.intron_ends);
            std::vector<int32_t>().swap(g.read_starts);
            std::vector<int32_t>().swap(g.read_ends);
        }

        // Sort Spliced RCs by (seq_id, rc_start, rc_end), with exact dplyr group_by tie breaking
        std::sort(sp_groups.begin(), sp_groups.end(), [](const SplicedRcGroup &a, const SplicedRcGroup &b) {
            if (a.key.seq_id != b.key.seq_id) return a.key.seq_id < b.key.seq_id;
            if (a.rc_start != b.rc_start) return a.rc_start < b.rc_start;
            if (a.rc_end != b.rc_end) return a.rc_end < b.rc_end;
            if (a.key.strand != b.key.strand) {
                return strand_factor_rank(a.key.strand) < strand_factor_rank(b.key.strand);
            }
            if (a.intron_ends_str != b.intron_ends_str) return a.intron_ends_str < b.intron_ends_str;
            if (a.intron_starts_str != b.intron_starts_str) return a.intron_starts_str < b.intron_starts_str;
            return a.key.conf_code < b.key.conf_code;
        });
    }

    // -------------------------------------------------------------------------
    // 2. Construct Unspliced ReadClasses
    // -------------------------------------------------------------------------
    std::vector<UnsplicedGroupRec> un_within_groups, un_new_groups;
    if (!bsd->unspliced.empty()) {
        // 2a. Build Reference Exons (from annotations + HC non-* spliced RCs)
        std::vector<RefExon> ref_exons;
        if (ann && ann->n_tx > 0) {
            std::unordered_map<std::string, int32_t> bsd_seq_map;
            for (size_t i = 0; i < bsd->seq_names.size(); ++i) {
                bsd_seq_map.emplace(bsd->seq_names[i], static_cast<int32_t>(i));
            }
            std::vector<int32_t> ann_to_bsd(ann->n_seq, -1);
            for (int32_t i = 0; i < ann->n_seq; ++i) {
                auto it = bsd_seq_map.find(ann->seq_names[i]);
                if (it != bsd_seq_map.end()) ann_to_bsd[i] = it->second;
            }
            ref_exons.reserve(ann->exon_start.size() + sp_groups.size() * 2);
            for (int32_t tx = 0; tx < ann->n_tx; ++tx) {
                int32_t sid = ann_to_bsd[ann->tx_seqid[tx]];
                if (sid < 0) continue;
                int8_t st = ann->tx_strand[tx];
                int32_t b = ann->tx_exon_ptr[tx];
                int32_t e = ann->tx_exon_ptr[tx + 1];
                for (int32_t k = b; k < e; ++k) {
                    ref_exons.push_back(RefExon{sid, ann->exon_start[k], ann->exon_end[k], st});
                }
            }
        }
        for (const auto &g : sp_groups) {
            if (g.key.conf_code == 0 && g.key.strand != 0) {
                size_t m = g.key.intron_starts.size();
                ref_exons.push_back(RefExon{g.key.seq_id, g.rc_start, g.key.intron_starts[0] - 1, g.key.strand});
                for (size_t k = 1; k < m; ++k) {
                    ref_exons.push_back(RefExon{g.key.seq_id, g.key.intron_ends[k - 1] + 1, g.key.intron_starts[k] - 1, g.key.strand});
                }
                ref_exons.push_back(RefExon{g.key.seq_id, g.key.intron_ends[m - 1] + 1, g.rc_end, g.key.strand});
            }
        }
        std::sort(ref_exons.begin(), ref_exons.end(), [](const RefExon &a, const RefExon &b) {
            if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
            if (a.start != b.start) return a.start < b.start;
            if (a.end != b.end) return a.end < b.end;
            return a.strand < b.strand;
        });
        ref_exons.erase(std::unique(ref_exons.begin(), ref_exons.end(), [](const RefExon &a, const RefExon &b) {
            return a.seq_id == b.seq_id && a.start == b.start && a.end == b.end && a.strand == b.strand;
        }), ref_exons.end());

        // 2b. Deduplicate single-exon reads by (seq_id, start, end, pos_strand) preserving encounter order
        struct UnspKey {
            int32_t seq_id, start, end;
            int8_t pos_strand;
            bool operator==(const UnspKey &o) const {
                return seq_id == o.seq_id && start == o.start && end == o.end && pos_strand == o.pos_strand;
            }
        };
        struct UnspKeyHash {
            size_t operator()(const UnspKey &k) const {
                uint64_t x = static_cast<uint64_t>(static_cast<uint32_t>(k.seq_id));
                x = splitmix64(x ^ static_cast<uint32_t>(k.start));
                x = splitmix64(x ^ static_cast<uint32_t>(k.end));
                return static_cast<size_t>(x ^ static_cast<uint64_t>(k.pos_strand));
            }
        };
        std::unordered_map<UnspKey, int32_t, UnspKeyHash> u_map;
        u_map.reserve(bsd->unspliced.size() / 2 + 64);
        std::vector<UniqueUnspliced> u_reads;
        u_reads.reserve(bsd->unspliced.size() / 2 + 64);
        for (const auto &ur : bsd->unspliced) {
            UnspKey k{ur.seq_id, ur.start, ur.end, ur.pos_strand};
            auto it = u_map.find(k);
            if (it == u_map.end()) {
                int32_t idx = static_cast<int32_t>(u_reads.size());
                u_map.emplace(k, idx);
                UniqueUnspliced rec{ur.seq_id, ur.start, ur.end, static_cast<int8_t>(ur.pos_strand ? 1 : -1), 1, {ur.read_id}};
                u_reads.push_back(std::move(rec));
            } else {
                auto &rec = u_reads[it->second];
                rec.count++;
                rec.read_ids.push_back(ur.read_id);
            }
        }
        std::unordered_map<UnspKey, int32_t, UnspKeyHash>().swap(u_map);

        std::vector<int32_t> all_u_idx(u_reads.size());
        std::iota(all_u_idx.begin(), all_u_idx.end(), 0);
        std::vector<uint8_t> matched_within(u_reads.size(), 0);

        // First pass: unsplicedWithin against referenceExons
        build_unspliced_classes_by_ref(
            u_reads, all_u_idx, ref_exons, stranded, "unsplicedWithin", matched_within, un_within_groups);

        // Second pass: remaining non-ref unspliced -> reduce -> unsplicedNew
        std::vector<int32_t> nonref_u_idx;
        std::vector<RefExon> to_reduce;
        nonref_u_idx.reserve(u_reads.size());
        to_reduce.reserve(u_reads.size());
        for (size_t i = 0; i < u_reads.size(); ++i) {
            if (!matched_within[i]) {
                nonref_u_idx.push_back(static_cast<int32_t>(i));
                const auto &u = u_reads[i];
                to_reduce.push_back(RefExon{u.seq_id, u.start, u.end, stranded ? u.strand : static_cast<int8_t>(0)});
            }
        }

        if (!nonref_u_idx.empty()) {
            std::sort(to_reduce.begin(), to_reduce.end(), [](const RefExon &a, const RefExon &b) {
                if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
                if (a.strand != b.strand) return strand_factor_rank(a.strand) < strand_factor_rank(b.strand);
                if (a.start != b.start) return a.start < b.start;
                return a.end < b.end;
            });
            std::vector<RefExon> reduced;
            reduced.reserve(to_reduce.size() / 2 + 16);
            for (const auto &cur : to_reduce) {
                if (reduced.empty() ||
                    reduced.back().seq_id != cur.seq_id ||
                    reduced.back().strand != cur.strand ||
                    cur.start > reduced.back().end + 1) {
                    reduced.push_back(cur);
                } else if (cur.end > reduced.back().end) {
                    reduced.back().end = cur.end;
                }
            }
            // Sort reduced by (seq_id, start, end) for plane sweep!
            std::sort(reduced.begin(), reduced.end(), [](const RefExon &a, const RefExon &b) {
                if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
                if (a.start != b.start) return a.start < b.start;
                return a.end < b.end;
            });
            std::vector<uint8_t> matched_new(u_reads.size(), 0);
            build_unspliced_classes_by_ref(
                u_reads, nonref_u_idx, reduced, stranded, "unsplicedNew", matched_new, un_new_groups);
        }
    }

    // -------------------------------------------------------------------------
    // 3. Serialize combined ReadClass list (Spliced + UnsplicedWithin + UnsplicedNew) to R
    // -------------------------------------------------------------------------
    int32_t n_sp_rc = static_cast<int32_t>(sp_groups.size());
    int32_t n_uw_rc = static_cast<int32_t>(un_within_groups.size());
    int32_t n_un_rc = static_cast<int32_t>(un_new_groups.size());
    int32_t n_tot_rc = n_sp_rc + n_uw_rc + n_un_rc;

    size_t tot_ex = 0;
    size_t tot_rid = 0;
    for (const auto &g : sp_groups) {
        tot_ex += g.key.intron_starts.size() + 1;
        tot_rid += g.read_ids.size();
    }
    tot_ex += static_cast<size_t>(n_uw_rc + n_un_rc);
    for (const auto &g : un_within_groups) tot_rid += g.read_ids.size();
    for (const auto &g : un_new_groups) tot_rid += g.read_ids.size();

    Rcpp::CharacterVector rc_names(n_tot_rc);
    Rcpp::IntegerVector chr_rc_idx(n_tot_rc);
    Rcpp::IntegerVector strand_code(n_tot_rc);
    Rcpp::NumericVector start_sd(n_tot_rc);
    Rcpp::NumericVector end_sd(n_tot_rc);
    Rcpp::IntegerVector read_count(n_tot_rc);
    Rcpp::IntegerVector read_count_pos(n_tot_rc);
    Rcpp::CharacterVector intron_starts(n_tot_rc);
    Rcpp::CharacterVector intron_ends(n_tot_rc);
    Rcpp::IntegerVector conf_code(n_tot_rc); // 1=HC, 2=LC, 3=unsplicedWithin, 4=unsplicedNew
    Rcpp::IntegerVector rc_ex_counts(n_tot_rc);
    Rcpp::IntegerVector rc_rid_counts(n_tot_rc);

    Rcpp::IntegerVector ex_seq_idx(tot_ex);
    Rcpp::IntegerVector ex_start(tot_ex);
    Rcpp::IntegerVector ex_end(tot_ex);
    Rcpp::IntegerVector ex_strand_code(tot_ex);
    Rcpp::IntegerVector ex_rank(tot_ex);
    Rcpp::IntegerVector ex_end_rank(tot_ex);
    Rcpp::IntegerVector flat_read_ids(tot_rid);

    size_t ex_pos = 0;
    size_t rid_pos = 0;
    int32_t rc_idx = 0;

    // Emit Spliced
    for (int32_t i = 0; i < n_sp_rc; ++i, ++rc_idx) {
        const auto &g = sp_groups[i];
        rc_names[rc_idx] = "rc." + std::to_string(i + 1);
        chr_rc_idx[rc_idx] = g.key.seq_id + 1;
        strand_code[rc_idx] = g.key.strand;
        start_sd[rc_idx] = g.start_sd;
        end_sd[rc_idx] = g.end_sd;
        read_count[rc_idx] = g.read_count;
        read_count_pos[rc_idx] = g.pos_strand_cnt;
        intron_starts[rc_idx] = g.intron_starts_str;
        intron_ends[rc_idx] = g.intron_ends_str;
        conf_code[rc_idx] = (g.key.conf_code == 0) ? 1 : 2;

        int32_t m = static_cast<int32_t>(g.key.intron_starts.size());
        int32_t n_ex = m + 1;
        rc_ex_counts[rc_idx] = n_ex;
        int32_t seq1 = g.key.seq_id + 1;
        int8_t st = g.key.strand;
        for (int32_t k = 0; k < n_ex; ++k) {
            int32_t es = (k == 0) ? g.rc_start : (g.key.intron_ends[k - 1] + 1);
            int32_t ee = (k == m) ? g.rc_end : (g.key.intron_starts[k] - 1);
            ex_seq_idx[ex_pos] = seq1;
            ex_start[ex_pos] = es;
            ex_end[ex_pos] = ee;
            ex_strand_code[ex_pos] = st;
            if (st == -1) {
                ex_rank[ex_pos] = n_ex - k;
                ex_end_rank[ex_pos] = k + 1;
            } else {
                ex_rank[ex_pos] = k + 1;
                ex_end_rank[ex_pos] = n_ex - k;
            }
            ex_pos++;
        }

        rc_rid_counts[rc_idx] = static_cast<int32_t>(g.read_ids.size());
        std::memcpy(flat_read_ids.begin() + rid_pos, g.read_ids.data(), g.read_ids.size() * sizeof(int32_t));
        rid_pos += g.read_ids.size();
    }

    auto emit_unspliced_vec = [&](const std::vector<UnsplicedGroupRec> &vec, int ccode) {
        for (const auto &g : vec) {
            rc_names[rc_idx] = g.rc_name;
            chr_rc_idx[rc_idx] = g.seq_id + 1;
            strand_code[rc_idx] = g.strand;
            start_sd[rc_idx] = g.start_sd;
            end_sd[rc_idx] = g.end_sd;
            read_count[rc_idx] = g.read_count;
            read_count_pos[rc_idx] = g.pos_strand_cnt;
            SET_STRING_ELT(intron_starts, rc_idx, NA_STRING);
            SET_STRING_ELT(intron_ends, rc_idx, NA_STRING);
            conf_code[rc_idx] = ccode;
            rc_ex_counts[rc_idx] = 1;

            ex_seq_idx[ex_pos] = g.seq_id + 1;
            ex_start[ex_pos] = g.start;
            ex_end[ex_pos] = g.end;
            ex_strand_code[ex_pos] = g.strand;
            ex_rank[ex_pos] = 1;
            ex_end_rank[ex_pos] = 1;
            ex_pos++;

            rc_rid_counts[rc_idx] = static_cast<int32_t>(g.read_ids.size());
            std::memcpy(flat_read_ids.begin() + rid_pos, g.read_ids.data(), g.read_ids.size() * sizeof(int32_t));
            rid_pos += g.read_ids.size();
            rc_idx++;
        }
    };
    emit_unspliced_vec(un_within_groups, 3);
    emit_unspliced_vec(un_new_groups, 4);

    // Optional tracked read names
    Rcpp::List tracked_info = R_NilValue;
    if (bsd->track_reads) {
        Rcpp::CharacterVector tnames(bsd->tracked_read_names.size());
        Rcpp::IntegerVector tids(bsd->tracked_read_ids.size());
        for (size_t i = 0; i < bsd->tracked_read_names.size(); ++i) {
            tnames[i] = bsd->tracked_read_names[i];
            tids[i] = bsd->tracked_read_ids[i];
        }
        tracked_info = Rcpp::List::create(Rcpp::Named("readNames") = tnames, Rcpp::Named("readId") = tids);
    }

    return Rcpp::List::create(
        Rcpp::Named("seq_names") = Rcpp::wrap(bsd->seq_names),
        Rcpp::Named("ref_seq_names") = Rcpp::wrap(bsd->ref_seq_names),
        Rcpp::Named("all_bam_sq_in_ref") = bsd->all_bam_sq_in_ref,
        Rcpp::Named("warnings") = Rcpp::wrap(bsd->warnings),
        Rcpp::Named("tracked") = tracked_info,
        Rcpp::Named("rc_meta") = Rcpp::List::create(
            Rcpp::Named("rc_names") = rc_names,
            Rcpp::Named("chr_rc_idx") = chr_rc_idx,
            Rcpp::Named("strand_code") = strand_code,
            Rcpp::Named("startSD") = start_sd,
            Rcpp::Named("endSD") = end_sd,
            Rcpp::Named("readCount") = read_count,
            Rcpp::Named("readCount_posStrand") = read_count_pos,
            Rcpp::Named("intronStarts") = intron_starts,
            Rcpp::Named("intronEnds") = intron_ends,
            Rcpp::Named("conf_code") = conf_code
        ),
        Rcpp::Named("exons") = Rcpp::List::create(
            Rcpp::Named("rc_ex_counts") = rc_ex_counts,
            Rcpp::Named("ex_seq_idx") = ex_seq_idx,
            Rcpp::Named("ex_start") = ex_start,
            Rcpp::Named("ex_end") = ex_end,
            Rcpp::Named("ex_strand_code") = ex_strand_code,
            Rcpp::Named("ex_rank") = ex_rank,
            Rcpp::Named("ex_end_rank") = ex_end_rank
        ),
        Rcpp::Named("reads") = Rcpp::List::create(
            Rcpp::Named("rc_rid_counts") = rc_rid_counts,
            Rcpp::Named("flat_read_ids") = flat_read_ids
        )
    );
    END_RCPP
}

extern "C" SEXP bambu_next_get_bsd_active_seqs_cpp(SEXP bsd_ptr_sexp) {
    BEGIN_RCPP
    Rcpp::XPtr<bambu_next::BamSampleData> bsd(bsd_ptr_sexp);
    size_t n_seq = bsd->seq_names.size();
    std::vector<uint8_t> active(n_seq, 0);
    for (int32_t r = 0; r < bsd->spliced.n_reads; ++r) {
        int32_t s = bsd->spliced.seq_id[r];
        if (s >= 0 && static_cast<size_t>(s) < n_seq) active[s] = 1;
    }
    for (const auto &ur : bsd->unspliced) {
        if (ur.seq_id >= 0 && static_cast<size_t>(ur.seq_id) < n_seq) active[ur.seq_id] = 1;
    }
    std::vector<int32_t> out;
    out.reserve(n_seq);
    for (size_t s = 0; s < n_seq; ++s) {
        if (active[s]) out.push_back(static_cast<int32_t>(s + 1));
    }
    return Rcpp::wrap(out);
    END_RCPP
}

extern "C" SEXP bambu_next_subset_bsd_seq_cpp(SEXP bsd_ptr_sexp, SEXP seq_1b_sexp) {
    BEGIN_RCPP
    Rcpp::XPtr<bambu_next::BamSampleData> bsd(bsd_ptr_sexp);
    int32_t target_sid = Rcpp::as<int32_t>(seq_1b_sexp) - 1;
    auto *sub = new bambu_next::BamSampleData();
    sub->run_name = bsd->run_name;
    sub->is_indexed = bsd->is_indexed;
    sub->track_reads = false;
    sub->all_bam_sq_in_ref = bsd->all_bam_sq_in_ref;
    sub->seq_names = bsd->seq_names;
    sub->ref_seq_names = bsd->ref_seq_names;
    sub->genome_seq_lengths = bsd->genome_seq_lengths;

    size_t n_j = bsd->unique_junctions.size();
    std::vector<int32_t> j_remap(n_j, -1);
    for (size_t j = 0; j < n_j; ++j) {
        if (bsd->unique_junctions[j].seq_id == target_sid) {
            j_remap[j] = static_cast<int32_t>(sub->unique_junctions.size());
            UniqueJunctionRec rec = bsd->unique_junctions[j];
            rec.merged_hc_id = -1;
            rec.merged_hc_id_no_na = -1;
            sub->unique_junctions.push_back(rec);
        }
    }

    int32_t n_sp = bsd->spliced.n_reads;
    sub->spliced.junc_ptr.push_back(0);
    for (int32_t r = 0; r < n_sp; ++r) {
        if (bsd->spliced.seq_id[r] != target_sid) continue;
        sub->spliced.seq_id.push_back(target_sid);
        sub->spliced.start.push_back(bsd->spliced.start[r]);
        sub->spliced.end.push_back(bsd->spliced.end[r]);
        sub->spliced.read_id.push_back(bsd->spliced.read_id[r]);
        sub->spliced.pos_strand.push_back(bsd->spliced.pos_strand[r]);
        int32_t jb = bsd->spliced.junc_ptr[r];
        int32_t je = bsd->spliced.junc_ptr[r + 1];
        for (int32_t k = jb; k < je; ++k) {
            int32_t nj = j_remap[bsd->spliced.junc_ids[k]];
            if (nj >= 0) sub->spliced.junc_ids.push_back(nj);
        }
        sub->spliced.junc_ptr.push_back(static_cast<int32_t>(sub->spliced.junc_ids.size()));
        sub->spliced.n_reads++;
    }

    for (const auto &ur : bsd->unspliced) {
        if (ur.seq_id == target_sid) sub->unspliced.push_back(ur);
    }
    return Rcpp::XPtr<bambu_next::BamSampleData>(sub, true);
    END_RCPP
}

extern "C" SEXP bambu_next_get_bsd_tracked_cpp(SEXP bsd_ptr_sexp) {
    BEGIN_RCPP
    Rcpp::XPtr<bambu_next::BamSampleData> bsd(bsd_ptr_sexp);
    if (!bsd->track_reads) return R_NilValue;
    Rcpp::CharacterVector tnames(bsd->tracked_read_names.size());
    Rcpp::IntegerVector tids(bsd->tracked_read_ids.size());
    for (size_t i = 0; i < bsd->tracked_read_names.size(); ++i) {
        tnames[i] = bsd->tracked_read_names[i];
        tids[i] = bsd->tracked_read_ids[i];
    }
    return Rcpp::List::create(Rcpp::Named("readNames") = tnames, Rcpp::Named("readId") = tids);
    END_RCPP
}

} // namespace bambu_next
