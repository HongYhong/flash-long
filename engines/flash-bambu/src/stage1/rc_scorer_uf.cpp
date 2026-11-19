#include <Rcpp.h>
#include <htslib/faidx.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/flat_csr_index.hpp"

namespace bambu_next {
namespace {

static inline int strand_rank_rc(int8_t s) {
    if (s == 1) return 1;
    if (s == -1) return 2;
    return 3;
}

static inline int strand_to_bucket(int8_t s) {
    if (s == 1) return 0;
    if (s == -1) return 1;
    return 2;
}

struct ExonIntervalSt {
    int32_t seq_id;
    int32_t start;
    int32_t end;
    int8_t strand;
};

struct GeneExonIndex {
    std::vector<std::string> gene_names; // sorted lexicographically ascending
    std::vector<int32_t> gene_ex_ptr;    // size n_gene + 1
    std::vector<ExonIntervalSt> exons;   // per gene sorted by (seq_id, start, end, strand) & reduced
    std::vector<int32_t> gene_seqid;     // primary seqid (or -1 if multi-chr? check per exon!)
    // Spatial index per (seq_id, strand_bucket): each gene interval -> gene_idx
    struct SpEntry {
        int32_t start;
        int32_t end;
        int32_t gene_idx;
    };
    struct Bucket {
        static constexpr int32_t BLOCK = 64;
        std::vector<SpEntry> entries;
        std::vector<int32_t> block_max_end;
        void finalize() {
            // Canonical IRanges NCList traversal order: (start asc, end desc, gene_idx asc)
            std::sort(entries.begin(), entries.end(), [](const SpEntry &a, const SpEntry &b) {
                if (a.start != b.start) return a.start < b.start;
                if (a.end != b.end) return a.end > b.end;
                return a.gene_idx < b.gene_idx;
            });
            size_t nb = (entries.size() + BLOCK - 1) / BLOCK;
            block_max_end.assign(nb, INT32_MIN);
            for (size_t i = 0; i < entries.size(); ++i) {
                size_t b = i / BLOCK;
                if (entries[i].end > block_max_end[b]) block_max_end[b] = entries[i].end;
            }
        }
        template <typename Fn>
        void query(int32_t q_start, int32_t q_end, int32_t min_ov, Fn &&cb) const {
            if (entries.empty() || (q_end - q_start + 1) < min_ov) return;
            // Want min(e.end, q_end) - max(e.start, q_start) + 1 >= min_ov
            // => e.start <= q_end - min_ov + 1 AND e.end >= q_start + min_ov - 1
            int32_t max_e_start = q_end - min_ov + 1;
            int32_t min_e_end = q_start + min_ov - 1;
            auto it = std::upper_bound(entries.begin(), entries.end(), max_e_start,
                [](int32_t val, const SpEntry &e) { return val < e.start; });
            int32_t limit = static_cast<int32_t>(it - entries.begin());
            if (limit <= 0) return;
            int32_t last_b = (limit - 1) / BLOCK;
            for (int32_t b = 0; b <= last_b; ++b) {
                if (block_max_end[b] < min_e_end) continue;
                int32_t beg = b * BLOCK;
                int32_t end = std::min(limit, beg + BLOCK);
                for (int32_t i = beg; i < end; ++i) {
                    if (entries[i].end >= min_e_end) {
                        int32_t ov = std::min(entries[i].end, q_end) - std::max(entries[i].start, q_start) + 1;
                        if (ov >= min_ov) cb(entries[i]);
                    }
                }
            }
        }
    };
    int32_t n_seq = 0;
    std::vector<Bucket> buckets; // size n_seq * 3
};

static void build_gene_exon_index_from_flat(
    int32_t n_seq,
    const std::vector<std::string> &tx_gene_names,
    const std::vector<int32_t> &tx_ex_counts,
    const std::vector<int32_t> &ex_seq_0b,
    const std::vector<int32_t> &ex_start,
    const std::vector<int32_t> &ex_end,
    const std::vector<int8_t> &ex_strand,
    GeneExonIndex &out
) {
    out = GeneExonIndex{};
    out.n_seq = n_seq;
    size_t n_tx = tx_gene_names.size();
    if (n_tx == 0) return;

    std::unordered_map<std::string, int32_t> gmap;
    gmap.reserve(n_tx / 2 + 16);
    std::vector<std::string> raw_genes;
    std::vector<int32_t> tx_raw_g(n_tx);
    for (size_t i = 0; i < n_tx; ++i) {
        auto it = gmap.find(tx_gene_names[i]);
        if (it == gmap.end()) {
            int32_t id = static_cast<int32_t>(raw_genes.size());
            gmap.emplace(tx_gene_names[i], id);
            raw_genes.push_back(tx_gene_names[i]);
            tx_raw_g[i] = id;
        } else {
            tx_raw_g[i] = it->second;
        }
    }

    int32_t n_gene = static_cast<int32_t>(raw_genes.size());
    std::vector<int32_t> perm(n_gene);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](int32_t a, int32_t b) {
        return raw_genes[a] < raw_genes[b];
    });
    std::vector<int32_t> raw_to_sorted(n_gene);
    out.gene_names.resize(n_gene);
    for (int32_t i = 0; i < n_gene; ++i) {
        raw_to_sorted[perm[i]] = i;
        out.gene_names[i] = std::move(raw_genes[perm[i]]);
    }

    std::vector<std::vector<ExonIntervalSt>> per_gene(n_gene);
    size_t pos = 0;
    for (size_t i = 0; i < n_tx; ++i) {
        int32_t gid = raw_to_sorted[tx_raw_g[i]];
        int32_t cnt = tx_ex_counts[i];
        for (int32_t k = 0; k < cnt; ++k, ++pos) {
            per_gene[gid].push_back(ExonIntervalSt{ex_seq_0b[pos], ex_start[pos], ex_end[pos], ex_strand[pos]});
        }
    }

    out.gene_ex_ptr.resize(n_gene + 1);
    out.buckets.resize(n_seq * 3);
    for (int32_t g = 0; g < n_gene; ++g) {
        out.gene_ex_ptr[g] = static_cast<int32_t>(out.exons.size());
        auto &v = per_gene[g];
        std::sort(v.begin(), v.end(), [](const ExonIntervalSt &a, const ExonIntervalSt &b) {
            if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
            if (a.strand != b.strand) return strand_rank_rc(a.strand) < strand_rank_rc(b.strand);
            if (a.start != b.start) return a.start < b.start;
            return a.end < b.end;
        });
        for (const auto &cur : v) {
            if (out.exons.size() == static_cast<size_t>(out.gene_ex_ptr[g]) ||
                out.exons.back().seq_id != cur.seq_id ||
                out.exons.back().strand != cur.strand ||
                cur.start > out.exons.back().end + 1) {
                out.exons.push_back(cur);
                if (cur.seq_id >= 0 && cur.seq_id < n_seq) {
                    out.buckets[cur.seq_id * 3 + strand_to_bucket(cur.strand)].entries.push_back(
                        GeneExonIndex::SpEntry{cur.start, cur.end, g});
                }
            } else if (cur.end > out.exons.back().end) {
                out.exons.back().end = cur.end;
                if (cur.seq_id >= 0 && cur.seq_id < n_seq) {
                    out.buckets[cur.seq_id * 3 + strand_to_bucket(cur.strand)].entries.back().end = cur.end;
                }
            }
        }
        std::vector<ExonIntervalSt>().swap(v);
    }
    out.gene_ex_ptr[n_gene] = static_cast<int32_t>(out.exons.size());
    for (auto &b : out.buckets) b.finalize();
}

static void build_gene_exon_index_from_csr(
    const FlatCSRAnnotation &csr,
    const std::unordered_map<std::string, int32_t> &seq_map,
    int32_t n_seq,
    GeneExonIndex &out
) {
    std::vector<int32_t> tx_ex_cnts(csr.n_tx);
    std::vector<std::string> tx_gnames(csr.n_tx);
    std::vector<int32_t> ex_seq_0b(csr.exon_start.size());
    std::vector<int8_t> ex_st(csr.exon_start.size());
    std::vector<int32_t> csr_seq_to_q(csr.n_seq, -1);
    for (int32_t i = 0; i < csr.n_seq; ++i) {
        auto it = seq_map.find(csr.seq_names[i]);
        if (it != seq_map.end()) csr_seq_to_q[i] = it->second;
    }
    for (int32_t tx = 0; tx < csr.n_tx; ++tx) {
        int32_t b = csr.tx_exon_ptr[tx];
        int32_t e = csr.tx_exon_ptr[tx + 1];
        tx_ex_cnts[tx] = e - b;
        tx_gnames[tx] = csr.gene_names[csr.tx_geneid[tx]];
        int32_t qsid = csr_seq_to_q[csr.tx_seqid[tx]];
        int8_t st = csr.tx_strand[tx];
        for (int32_t k = b; k < e; ++k) {
            ex_seq_0b[k] = qsid;
            ex_st[k] = st;
        }
    }
    build_gene_exon_index_from_flat(
        n_seq, tx_gnames, tx_ex_cnts, ex_seq_0b, csr.exon_start, csr.exon_end, ex_st, out);
}

// Compute IRanges::intersect total width between query's exons and gene's reduced exons
static int64_t compute_intersect_width(
    const std::vector<std::pair<int32_t, int32_t>> &q_reduced_ranges,
    const GeneExonIndex &idx,
    int32_t gene_idx
) {
    // Gene's intervals across all seqids/strands vs query's ranges (IRanges::intersect strips strand & chr?
    // Wait: lines 449-450 of bambu use intersect(ranges(grl), ranges(geneRanges)), which operates on pure IRanges!)
    // Since ov already filtered chr & strand compatibility, normally all exons of gene are on same chr/strand,
    // EXCEPT if gene has mixed strand or multi-seq (rare). Replicating exact IRanges::intersect(ranges(q), ranges(gene)):
    int32_t b = idx.gene_ex_ptr[gene_idx];
    int32_t e = idx.gene_ex_ptr[gene_idx + 1];
    std::vector<std::pair<int32_t, int32_t>> g_rng;
    g_rng.reserve(e - b);
    for (int32_t k = b; k < e; ++k) g_rng.emplace_back(idx.exons[k].start, idx.exons[k].end);
    std::sort(g_rng.begin(), g_rng.end());
    // Merge g_rng into disjoint sorted ranges
    size_t w = 0;
    for (size_t i = 0; i < g_rng.size(); ++i) {
        if (w == 0 || g_rng[i].first > g_rng[w - 1].second + 1) {
            g_rng[w++] = g_rng[i];
        } else if (g_rng[i].second > g_rng[w - 1].second) {
            g_rng[w - 1].second = g_rng[i].second;
        }
    }
    g_rng.resize(w);

    int64_t tot = 0;
    size_t i = 0, j = 0;
    while (i < q_reduced_ranges.size() && j < g_rng.size()) {
        int32_t s = std::max(q_reduced_ranges[i].first, g_rng[j].first);
        int32_t en = std::min(q_reduced_ranges[i].second, g_rng[j].second);
        if (s <= en) tot += static_cast<int64_t>(en - s + 1);
        if (q_reduced_ranges[i].second < g_rng[j].second) i++;
        else j++;
    }
    return tot;
}

static void assign_gene_ids_by_ref(
    const std::vector<int32_t> &q_subset,
    const std::vector<int32_t> &q_ex_ptr,
    const std::vector<int32_t> &q_ex_seq,
    const std::vector<int32_t> &q_ex_start,
    const std::vector<int32_t> &q_ex_end,
    const std::vector<int8_t> &q_ex_strand,
    const GeneExonIndex &ref_idx,
    int32_t min_overlap,
    bool fusion_mode,
    std::vector<std::string> &out_gene_ids
) {
    if (q_subset.empty() || ref_idx.gene_names.empty()) return;
    struct HitKey {
        int32_t gid;
        int32_t ov;
        int32_t q_ex_k;
        int32_t s_start;
        int32_t s_end;
    };
    std::vector<HitKey> hit_ov;
    struct CandGene {
        int32_t gid;
        int32_t first_k;
        int32_t first_s_start;
        int32_t first_s_end;
    };
    std::vector<CandGene> hit_genes;
    std::vector<std::pair<int32_t, int32_t>> q_rng;
    hit_ov.reserve(16);
    hit_genes.reserve(16);
    q_rng.reserve(32);

    auto rank_less = [](const CandGene &a, const CandGene &b) {
        if (a.first_k != b.first_k) return a.first_k < b.first_k;
        if (a.first_s_start != b.first_s_start) return a.first_s_start < b.first_s_start;
        if (a.first_s_end != b.first_s_end) return a.first_s_end > b.first_s_end;
        return a.gid < b.gid;
    };

    // Count per-chromosome query exons vs subject exons to match IRanges C_find_overlaps_in_groups_NCList
    // (which builds NCList on query when n_q_chr < n_s_chr -> ordering ties by gid asc,
    //  and builds NCList on subject when n_q_chr >= n_s_chr -> ordering ties by s_start asc, s_end desc, gid asc)
    std::vector<int64_t> n_q_chr(ref_idx.n_seq, 0);
    std::vector<int64_t> n_s_chr(ref_idx.n_seq, 0);
    for (int32_t s = 0; s < ref_idx.n_seq; ++s) {
        n_s_chr[s] = static_cast<int64_t>(ref_idx.buckets[s * 3 + 0].entries.size() +
                                          ref_idx.buckets[s * 3 + 1].entries.size() +
                                          ref_idx.buckets[s * 3 + 2].entries.size());
    }
    for (int32_t q : q_subset) {
        for (int32_t k = q_ex_ptr[q]; k < q_ex_ptr[q + 1]; ++k) {
            int32_t sid = q_ex_seq[k];
            if (sid >= 0 && sid < ref_idx.n_seq) n_q_chr[sid]++;
        }
    }

    for (int32_t q : q_subset) {
        int32_t qb = q_ex_ptr[q];
        int32_t qe = q_ex_ptr[q + 1];
        if (qb >= qe) continue;
        int32_t q_sid = q_ex_seq[qb];
        bool nclist_on_q = (q_sid >= 0 && q_sid < ref_idx.n_seq) && (n_q_chr[q_sid] < n_s_chr[q_sid]);
        hit_ov.clear();

        for (int32_t k = qb; k < qe; ++k) {
            int32_t sid = q_ex_seq[k];
            if (sid < 0 || sid >= ref_idx.n_seq) continue;
            int8_t qst = q_ex_strand[k];
            int32_t qs = q_ex_start[k];
            int32_t qen = q_ex_end[k];
            int32_t base = sid * 3;
            int32_t rk_k = k - qb;
            auto add_cb = [&](const GeneExonIndex::SpEntry &e) {
                int32_t ov = std::min(e.end, qen) - std::max(e.start, qs) + 1;
                hit_ov.push_back(HitKey{e.gene_idx, ov, rk_k, e.start, e.end});
            };
            if (qst == 0) {
                ref_idx.buckets[base + 0].query(qs, qen, 1, add_cb);
                ref_idx.buckets[base + 1].query(qs, qen, 1, add_cb);
                ref_idx.buckets[base + 2].query(qs, qen, 1, add_cb);
            } else if (qst == 1) {
                ref_idx.buckets[base + 0].query(qs, qen, 1, add_cb);
                ref_idx.buckets[base + 2].query(qs, qen, 1, add_cb);
            } else {
                ref_idx.buckets[base + 1].query(qs, qen, 1, add_cb);
                ref_idx.buckets[base + 2].query(qs, qen, 1, add_cb);
            }
        }
        if (hit_ov.empty()) continue;
        std::sort(hit_ov.begin(), hit_ov.end(), [](const HitKey &a, const HitKey &b) {
            if (a.gid != b.gid) return a.gid < b.gid;
            if (a.q_ex_k != b.q_ex_k) return a.q_ex_k < b.q_ex_k;
            if (a.s_start != b.s_start) return a.s_start < b.s_start;
            return a.s_end > b.s_end;
        });
        hit_genes.clear();
        for (size_t i = 0; i < hit_ov.size(); ) {
            int32_t gid = hit_ov[i].gid;
            int64_t ov_sum = 0;
            CandGene cg{gid, hit_ov[i].q_ex_k, hit_ov[i].s_start, hit_ov[i].s_end};
            while (i < hit_ov.size() && hit_ov[i].gid == gid) {
                ov_sum += hit_ov[i].ov;
                ++i;
            }
            if (ov_sum >= min_overlap) hit_genes.push_back(cg);
        }
        if (hit_genes.empty()) continue;

        if (hit_genes.size() == 1) {
            out_gene_ids[q] = ref_idx.gene_names[hit_genes[0].gid];
            continue;
        }

        // Multiple hit genes: compute exact IRanges::intersect width
        q_rng.clear();
        for (int32_t k = qb; k < qe; ++k) q_rng.emplace_back(q_ex_start[k], q_ex_end[k]);
        std::sort(q_rng.begin(), q_rng.end());
        size_t w = 0;
        for (size_t i = 0; i < q_rng.size(); ++i) {
            if (w == 0 || q_rng[i].first > q_rng[w - 1].second + 1) {
                q_rng[w++] = q_rng[i];
            } else if (q_rng[i].second > q_rng[w - 1].second) {
                q_rng[w - 1].second = q_rng[i].second;
            }
        }
        q_rng.resize(w);

        // Order hit_genes matching Bioconductor IRanges::find_overlaps_in_groups_NCList tie order
        if (nclist_on_q) {
            std::sort(hit_genes.begin(), hit_genes.end(), [](const CandGene &a, const CandGene &b) {
                if (a.first_k != b.first_k) return a.first_k < b.first_k;
                return a.gid < b.gid;
            });
        } else {
            std::sort(hit_genes.begin(), hit_genes.end(), rank_less);
        }

        if (fusion_mode) {
            std::string fused;
            for (const auto &cg : hit_genes) {
                int64_t iw = compute_intersect_width(q_rng, ref_idx, cg.gid);
                if (iw > min_overlap) {
                    if (!fused.empty()) fused.push_back(':');
                    fused += ref_idx.gene_names[cg.gid];
                }
            }
            if (!fused.empty()) out_gene_ids[q] = std::move(fused);
        } else {
            int64_t best_w = -1;
            int32_t best_gid = -1;
            for (const auto &cg : hit_genes) {
                int64_t iw = compute_intersect_width(q_rng, ref_idx, cg.gid);
                if (iw > best_w) {
                    best_w = iw;
                    best_gid = cg.gid;
                }
            }
            if (best_gid >= 0) out_gene_ids[q] = ref_idx.gene_names[best_gid];
        }
    }
}

struct UnionFindMin {
    std::vector<int32_t> parent;
    std::vector<int32_t> min_val;
    void init(int32_t n) {
        parent.resize(n);
        min_val.resize(n);
        for (int32_t i = 0; i < n; ++i) {
            parent[i] = i;
            min_val[i] = i + 1; // 1-based newExonId!
        }
    }
    int32_t find(int32_t x) {
        int32_t r = x;
        while (parent[r] != r) r = parent[r];
        while (parent[x] != r) {
            int32_t nxt = parent[x];
            parent[x] = r;
            x = nxt;
        }
        return r;
    }
    void unite(int32_t a, int32_t b) {
        int32_t ra = find(a);
        int32_t rb = find(b);
        if (ra == rb) return;
        if (ra > rb) std::swap(ra, rb);
        parent[rb] = ra;
        if (min_val[rb] < min_val[ra]) min_val[ra] = min_val[rb];
    }
};

// Replicate assignGeneIdsNoReference exactly using O(M log M) sweep + Union-Find
static void assign_gene_ids_no_ref(
    const std::vector<int32_t> &q_subset,
    const std::vector<int32_t> &q_ex_ptr,
    const std::vector<int32_t> &q_ex_seq,
    const std::vector<int32_t> &q_ex_start,
    const std::vector<int32_t> &q_ex_end,
    const std::vector<int8_t> &q_ex_strand,
    const std::string &prefix,
    std::vector<std::string> &out_gene_ids
) {
    if (q_subset.empty()) return;
    struct ExOwner {
        int32_t seq_id;
        int8_t strand;
        int32_t start;
        int32_t end;
        int32_t sub_pos; // 0 .. q_subset.size()-1
    };
    std::vector<ExOwner> all_ex;
    for (size_t si = 0; si < q_subset.size(); ++si) {
        int32_t q = q_subset[si];
        for (int32_t k = q_ex_ptr[q]; k < q_ex_ptr[q + 1]; ++k) {
            all_ex.push_back(ExOwner{q_ex_seq[k], q_ex_strand[k], q_ex_start[k], q_ex_end[k], static_cast<int32_t>(si)});
        }
    }
    if (all_ex.empty()) return;

    // Sort matches GenomicRanges::reduce(unlist(grl)): (seq_id, strand [+, -, *], start, end)
    std::sort(all_ex.begin(), all_ex.end(), [](const ExOwner &a, const ExOwner &b) {
        if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
        if (a.strand != b.strand) return strand_rank_rc(a.strand) < strand_rank_rc(b.strand);
        if (a.start != b.start) return a.start < b.start;
        return a.end < b.end;
    });

    // Sweep to merge intervals & build edges to sub_pos (and link all merged_exon_ids owned by same tx!)
    std::vector<int32_t> tx_first_merged_exon(q_subset.size(), -1);
    int32_t n_merged = 0;
    int32_t cur_seq = -1, cur_st = 99, cur_end = -1;
    std::vector<std::pair<int32_t, int32_t>> edges; // (merged_exon_id_a, merged_exon_id_b)

    for (const auto &e : all_ex) {
        bool new_seg = (n_merged == 0 || e.seq_id != cur_seq || e.strand != cur_st || e.start > cur_end + 1);
        int32_t me_id;
        if (new_seg) {
            me_id = n_merged++;
            cur_seq = e.seq_id;
            cur_st = e.strand;
            cur_end = e.end;
        } else {
            me_id = n_merged - 1;
            if (e.end > cur_end) cur_end = e.end;
        }
        int32_t prev_me = tx_first_merged_exon[e.sub_pos];
        if (prev_me < 0) {
            tx_first_merged_exon[e.sub_pos] = me_id;
        } else if (prev_me != me_id) {
            edges.emplace_back(prev_me, me_id);
        }
    }

    UnionFindMin uf;
    uf.init(n_merged);
    for (const auto &ed : edges) {
        uf.unite(ed.first, ed.second);
    }

    std::string full_prefix = prefix + "Gene";
    for (size_t si = 0; si < q_subset.size(); ++si) {
        int32_t q = q_subset[si];
        int32_t me = tx_first_merged_exon[si];
        if (me >= 0) {
            int32_t r = uf.find(me);
            out_gene_ids[q] = full_prefix + std::to_string(uf.min_val[r]);
        }
    }
}

} // namespace
} // namespace bambu_next

// R-callable C++ assignGeneIds
extern "C" SEXP bambu_next_assign_gene_ids_cpp(
    SEXP q_seq_names_sexp,
    SEXP q_tx_ex_counts_sexp,
    SEXP q_ex_seq_1b_sexp,
    SEXP q_ex_start_sexp,
    SEXP q_ex_end_sexp,
    SEXP q_ex_strand_sexp,
    SEXP csr_ptr_sexp,
    SEXP ann_ex_counts_sexp,
    SEXP ann_gene_names_sexp,
    SEXP ann_ex_seq_1b_sexp,
    SEXP ann_ex_start_sexp,
    SEXP ann_ex_end_sexp,
    SEXP ann_ex_strand_sexp,
    SEXP min_exon_overlap_sexp,
    SEXP fusion_mode_sexp,
    SEXP max_chain_iter_sexp
) {
    BEGIN_RCPP
    Rcpp::CharacterVector q_snames(q_seq_names_sexp);
    Rcpp::IntegerVector q_ex_cnts(q_tx_ex_counts_sexp);
    Rcpp::IntegerVector q_es(q_ex_seq_1b_sexp);
    Rcpp::IntegerVector q_s(q_ex_start_sexp);
    Rcpp::IntegerVector q_e(q_ex_end_sexp);
    Rcpp::IntegerVector q_st(q_ex_strand_sexp);

    int32_t min_ov = Rf_asInteger(min_exon_overlap_sexp);
    bool fusion_mode = Rcpp::as<bool>(fusion_mode_sexp);
    int32_t max_chain = Rf_asInteger(max_chain_iter_sexp);

    int32_t n_q = static_cast<int32_t>(q_ex_cnts.size());
    int32_t n_seq = static_cast<int32_t>(q_snames.size());
    std::unordered_map<std::string, int32_t> q_seq_map;
    q_seq_map.reserve(n_seq * 2 + 8);
    for (int32_t i = 0; i < n_seq; ++i) q_seq_map.emplace(Rcpp::as<std::string>(q_snames[i]), i);

    std::vector<int32_t> q_ex_ptr(n_q + 1);
    std::vector<int8_t> q_tx_st(n_q, 0);
    int32_t off = 0;
    for (int32_t i = 0; i < n_q; ++i) {
        q_ex_ptr[i] = off;
        off += q_ex_cnts[i];
        if (q_ex_cnts[i] > 0) q_tx_st[i] = static_cast<int8_t>(q_st[q_ex_ptr[i]]);
    }
    q_ex_ptr[n_q] = off;

    std::vector<int32_t> q_ex_seq_vec(off), q_ex_start_vec(off), q_ex_end_vec(off);
    std::vector<int8_t> q_ex_st_vec(off);
    for (int32_t i = 0; i < off; ++i) {
        q_ex_seq_vec[i] = q_es[i] - 1;
        q_ex_start_vec[i] = q_s[i];
        q_ex_end_vec[i] = q_e[i];
        q_ex_st_vec[i] = static_cast<int8_t>(q_st[i]);
    }

    // Build reference GeneExonIndex (from FlatCSRAnnotation if provided, or explicit flat annotation arrays)
    bambu_next::GeneExonIndex ref_idx;
    std::unordered_set<std::string> ref_gene_set;
    if (!Rf_isNull(csr_ptr_sexp)) {
        Rcpp::XPtr<bambu_next::FlatCSRAnnotation> csr(csr_ptr_sexp);
        if (csr->n_tx > 0) {
            bambu_next::build_gene_exon_index_from_csr(*csr, q_seq_map, n_seq, ref_idx);
            for (const auto &g : csr->gene_names) ref_gene_set.insert(g);
        }
    } else if (!Rf_isNull(ann_ex_counts_sexp)) {
        Rcpp::IntegerVector a_cnts(ann_ex_counts_sexp);
        if (a_cnts.size() > 0) {
            Rcpp::CharacterVector a_gn(ann_gene_names_sexp);
            Rcpp::IntegerVector a_es(ann_ex_seq_1b_sexp);
            Rcpp::IntegerVector a_s(ann_ex_start_sexp);
            Rcpp::IntegerVector a_e(ann_ex_end_sexp);
            Rcpp::IntegerVector a_st(ann_ex_strand_sexp);
            size_t na = a_cnts.size();
            std::vector<std::string> gn_v(na);
            std::vector<int32_t> cnt_v(na);
            for (size_t i = 0; i < na; ++i) {
                gn_v[i] = Rcpp::as<std::string>(a_gn[i]);
                ref_gene_set.insert(gn_v[i]);
                cnt_v[i] = a_cnts[i];
            }
            size_t nae = a_es.size();
            std::vector<int32_t> aes_v(nae), as_v(nae), ae_v(nae);
            std::vector<int8_t> ast_v(nae);
            for (size_t i = 0; i < nae; ++i) {
                aes_v[i] = a_es[i] - 1;
                as_v[i] = a_s[i];
                ae_v[i] = a_e[i];
                ast_v[i] = static_cast<int8_t>(a_st[i]);
            }
            bambu_next::build_gene_exon_index_from_flat(n_seq, gn_v, cnt_v, aes_v, as_v, ae_v, ast_v, ref_idx);
        }
    }

    std::vector<std::string> out_gene(n_q); // empty = NA
    std::vector<int32_t> stranded_idx, unstranded_idx;
    stranded_idx.reserve(n_q);
    for (int32_t i = 0; i < n_q; ++i) {
        if (q_tx_st[i] != 0) stranded_idx.push_back(i);
        else unstranded_idx.push_back(i);
    }

    bool has_colon_warning = false;
    if (!ref_idx.gene_names.empty()) {
        if (fusion_mode) {
            for (const auto &g : ref_idx.gene_names) {
                if (g.find(':') != std::string::npos) {
                    has_colon_warning = true;
                    break;
                }
            }
        }
        bambu_next::assign_gene_ids_by_ref(
            stranded_idx, q_ex_ptr, q_ex_seq_vec, q_ex_start_vec, q_ex_end_vec, q_ex_st_vec,
            ref_idx, min_ov, fusion_mode, out_gene);

        // Chaining iterations (up to max_chain passes total, matching lines 388-403)
        std::vector<int32_t> new_gene_set, ref_gene_subset;
        for (int32_t q : stranded_idx) {
            if (out_gene[q].empty()) new_gene_set.push_back(q);
            else ref_gene_subset.push_back(q);
        }
        auto build_sub_idx = [&](const std::vector<int32_t> &sub_q, bambu_next::GeneExonIndex &sub_idx) {
            std::vector<std::string> gn;
            std::vector<int32_t> cnts, esq, est, een;
            std::vector<int8_t> estr;
            gn.reserve(sub_q.size());
            cnts.reserve(sub_q.size());
            for (int32_t sq : sub_q) {
                gn.push_back(out_gene[sq]);
                int32_t b = q_ex_ptr[sq], e = q_ex_ptr[sq + 1];
                cnts.push_back(e - b);
                for (int32_t k = b; k < e; ++k) {
                    esq.push_back(q_ex_seq_vec[k]);
                    est.push_back(q_ex_start_vec[k]);
                    een.push_back(q_ex_end_vec[k]);
                    estr.push_back(q_ex_st_vec[k]);
                }
            }
            bambu_next::build_gene_exon_index_from_flat(n_seq, gn, cnts, esq, est, een, estr, sub_idx);
        };

        if (!new_gene_set.empty() && !ref_gene_subset.empty()) {
            bambu_next::GeneExonIndex sub_idx;
            build_sub_idx(ref_gene_subset, sub_idx);
            bambu_next::assign_gene_ids_by_ref(
                new_gene_set, q_ex_ptr, q_ex_seq_vec, q_ex_start_vec, q_ex_end_vec, q_ex_st_vec,
                sub_idx, min_ov, false, out_gene);

            int32_t chain_cnt = 1;
            while (chain_cnt < max_chain) {
                std::vector<int32_t> newly_assigned, still_unassigned;
                for (int32_t q : new_gene_set) {
                    if (!out_gene[q].empty()) newly_assigned.push_back(q);
                    else still_unassigned.push_back(q);
                }
                if (newly_assigned.empty() || still_unassigned.empty()) break;
                bambu_next::GeneExonIndex chain_idx;
                build_sub_idx(newly_assigned, chain_idx);
                bambu_next::assign_gene_ids_by_ref(
                    still_unassigned, q_ex_ptr, q_ex_seq_vec, q_ex_start_vec, q_ex_end_vec, q_ex_st_vec,
                    chain_idx, min_ov, false, out_gene);
                new_gene_set.swap(still_unassigned);
                chain_cnt++;
            }
        }
    }

    // Remaining NA among stranded -> BambuGene*
    {
        std::vector<int32_t> rem_stranded;
        for (int32_t q : stranded_idx) {
            if (out_gene[q].empty()) rem_stranded.push_back(q);
        }
        bambu_next::assign_gene_ids_no_ref(
            rem_stranded, q_ex_ptr, q_ex_seq_vec, q_ex_start_vec, q_ex_end_vec, q_ex_st_vec,
            "Bambu", out_gene);
    }

    // Unstranded queries -> match assigned RCs -> match annotations -> unstranded.Gene*
    if (!unstranded_idx.empty()) {
        std::vector<int32_t> assigned_q;
        for (int32_t i = 0; i < n_q; ++i) {
            if (!out_gene[i].empty()) assigned_q.push_back(i);
        }
        if (!assigned_q.empty()) {
            bambu_next::GeneExonIndex rc_idx;
            std::vector<std::string> gn;
            std::vector<int32_t> cnts, esq, est, een;
            std::vector<int8_t> estr;
            for (int32_t aq : assigned_q) {
                gn.push_back(out_gene[aq]);
                int32_t b = q_ex_ptr[aq], e = q_ex_ptr[aq + 1];
                cnts.push_back(e - b);
                for (int32_t k = b; k < e; ++k) {
                    esq.push_back(q_ex_seq_vec[k]);
                    est.push_back(q_ex_start_vec[k]);
                    een.push_back(q_ex_end_vec[k]);
                    estr.push_back(q_ex_st_vec[k]);
                }
            }
            bambu_next::build_gene_exon_index_from_flat(n_seq, gn, cnts, esq, est, een, estr, rc_idx);
            bambu_next::assign_gene_ids_by_ref(
                unstranded_idx, q_ex_ptr, q_ex_seq_vec, q_ex_start_vec, q_ex_end_vec, q_ex_st_vec,
                rc_idx, min_ov, false, out_gene);
        }
        if (!ref_idx.gene_names.empty()) {
            std::vector<int32_t> rem_u;
            for (int32_t q : unstranded_idx) {
                if (out_gene[q].empty()) rem_u.push_back(q);
            }
            bambu_next::assign_gene_ids_by_ref(
                rem_u, q_ex_ptr, q_ex_seq_vec, q_ex_start_vec, q_ex_end_vec, q_ex_st_vec,
                ref_idx, min_ov, false, out_gene);
        }
        std::vector<int32_t> rem_u2;
        for (int32_t q : unstranded_idx) {
            if (out_gene[q].empty()) rem_u2.push_back(q);
        }
        bambu_next::assign_gene_ids_no_ref(
            rem_u2, q_ex_ptr, q_ex_seq_vec, q_ex_start_vec, q_ex_end_vec, q_ex_st_vec,
            "unstranded.", out_gene);
    }

    Rcpp::CharacterVector res_gene(n_q);
    Rcpp::LogicalVector res_novel(n_q);
    for (int32_t i = 0; i < n_q; ++i) {
        if (out_gene[i].empty()) {
            SET_STRING_ELT(res_gene, i, NA_STRING);
            res_novel[i] = true;
        } else {
            res_gene[i] = out_gene[i];
            res_novel[i] = (ref_gene_set.find(out_gene[i]) == ref_gene_set.end());
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("GENEID") = res_gene,
        Rcpp::Named("novelGene") = res_novel,
        Rcpp::Named("colon_warning") = has_colon_warning
    );
    END_RCPP
}

// Fast C++ isReadClassCompatible against FlatCSRAnnotation
extern "C" SEXP bambu_next_is_rc_compatible_cpp(
    SEXP q_seq_names_sexp,
    SEXP q_tx_ex_counts_sexp,
    SEXP q_ex_seq_1b_sexp,
    SEXP q_ex_start_sexp,
    SEXP q_ex_end_sexp,
    SEXP q_ex_strand_sexp,
    SEXP q_ex_rank_sexp,
    SEXP q_ex_end_rank_sexp,
    SEXP csr_ptr_sexp,
    SEXP n_threads_sexp
) {
    BEGIN_RCPP
    Rcpp::IntegerVector q_ex_cnts(q_tx_ex_counts_sexp);
    int32_t n_q = static_cast<int32_t>(q_ex_cnts.size());
    Rcpp::IntegerVector out_comp(n_q, 0);
    Rcpp::LogicalVector out_equal(n_q, false);
    if (n_q == 0 || Rf_isNull(csr_ptr_sexp)) {
        return Rcpp::List::create(Rcpp::Named("compatible") = out_comp, Rcpp::Named("equal") = out_equal);
    }
    Rcpp::XPtr<bambu_next::FlatCSRAnnotation> csr(csr_ptr_sexp);
    if (csr->n_tx == 0) {
        return Rcpp::List::create(Rcpp::Named("compatible") = out_comp, Rcpp::Named("equal") = out_equal);
    }

    Rcpp::CharacterVector q_snames(q_seq_names_sexp);
    Rcpp::IntegerVector q_es(q_ex_seq_1b_sexp);
    Rcpp::IntegerVector q_s(q_ex_start_sexp);
    Rcpp::IntegerVector q_e(q_ex_end_sexp);
    Rcpp::IntegerVector q_st(q_ex_strand_sexp);
    Rcpp::IntegerVector q_rk(q_ex_rank_sexp);
    Rcpp::IntegerVector q_erk(q_ex_end_rank_sexp);
    int n_threads = std::max(1, Rcpp::as<int>(n_threads_sexp));

    // Map q seqnames -> csr seqid
    std::vector<int32_t> q_to_csr_seq(q_snames.size(), -1);
    for (R_xlen_t i = 0; i < q_snames.size(); ++i) {
        auto it = csr->seq_to_id.find(Rcpp::as<std::string>(q_snames[i]));
        if (it != csr->seq_to_id.end()) q_to_csr_seq[i] = it->second;
    }

    // Pre-build annotated intron lookup hash (per csr_seqid, start, end, strand)
    struct AnnIntKey {
        int32_t sid, s, e;
        int8_t st;
        bool operator==(const AnnIntKey &o) const {
            return sid == o.sid && s == o.s && e == o.e && st == o.st;
        }
    };
    struct AnnIntHash {
        size_t operator()(const AnnIntKey &k) const {
            uint64_t x = static_cast<uint64_t>(static_cast<uint32_t>(k.sid));
            x = (x * 0x9e3779b97f4a7c15ULL) ^ static_cast<uint32_t>(k.s);
            x = (x * 0xbf58476d1ce4e5b9ULL) ^ static_cast<uint32_t>(k.e);
            x = (x * 0x94d049bb133111ebULL) ^ static_cast<uint32_t>(k.st + 2);
            return static_cast<size_t>(x ^ (x >> 30));
        }
    };
    std::unordered_set<AnnIntKey, AnnIntHash> ann_introns;
    ann_introns.reserve(csr->exon_start.size() / 2 + 64);
    for (int32_t tx = 0; tx < csr->n_tx; ++tx) {
        int32_t b = csr->tx_exon_ptr[tx];
        int32_t e = csr->tx_exon_ptr[tx + 1];
        if (e - b < 2) continue;
        int32_t sid = csr->tx_seqid[tx];
        int8_t st = csr->tx_strand[tx];
        for (int32_t k = b; k + 1 < e; ++k) {
            int32_t ist = csr->exon_end[k] + 1;
            int32_t ien = csr->exon_start[k + 1] - 1;
            if (ist <= ien) ann_introns.insert(AnnIntKey{sid, ist, ien, st});
        }
    }

    std::vector<int32_t> q_ptr(n_q + 1);
    int32_t off = 0;
    for (int32_t i = 0; i < n_q; ++i) {
        q_ptr[i] = off;
        off += q_ex_cnts[i];
    }
    q_ptr[n_q] = off;

    // Copy flat cutStartEnd exons for queries into plain std::vectors so C++ threads never touch SEXPs!
    std::vector<int32_t> q_csr_sid(n_q, -1);
    std::vector<int8_t> q_strand(n_q, 0);
    std::vector<int32_t> q_cut_s(off), q_cut_e(off);
    std::vector<int32_t> q_cut_mins(n_q), q_cut_maxe(n_q);
    std::vector<uint8_t> q_pass_intron(n_q, 1);

    for (int32_t q = 0; q < n_q; ++q) {
        int32_t b = q_ptr[q];
        int32_t e = q_ptr[q + 1];
        if (b >= e) {
            q_pass_intron[q] = 0;
            continue;
        }
        int32_t csid = q_to_csr_seq[q_es[b] - 1];
        q_csr_sid[q] = csid;
        int8_t st = static_cast<int8_t>(q_st[b]);
        q_strand[q] = st;
        if (csid < 0) {
            q_pass_intron[q] = 0;
            continue;
        }
        // Intron prefilter on cutStartEndFromGrangesList(query)
        int32_t cmin = INT32_MAX, cmax = INT32_MIN;
        for (int32_t k = b; k < e; ++k) {
            q_cut_s[k] = q_s[k];
            q_cut_e[k] = q_e[k];
        }
        for (int32_t k = b; k < e; ++k) {
            bool is_start = (st != -1) ? (q_rk[k] == 1) : (q_erk[k] == 1);
            if (is_start) q_cut_s[k] = q_cut_e[k] - 1;
        }
        for (int32_t k = b; k < e; ++k) {
            bool is_end = (st == -1) ? (q_rk[k] == 1) : (q_erk[k] == 1);
            if (is_end) q_cut_e[k] = q_cut_s[k] + 1;
        }
        for (int32_t k = b; k < e; ++k) {
            if (q_cut_s[k] < cmin) cmin = q_cut_s[k];
            if (q_cut_e[k] > cmax) cmax = q_cut_e[k];
        }
        q_cut_mins[q] = cmin;
        q_cut_maxe[q] = cmax;

        // Check every gap in cut_query exists in ann_introns (exact strand equality!)
        for (int32_t k = b; k + 1 < e; ++k) {
            int32_t ist = q_cut_e[k] + 1;
            int32_t ien = q_cut_s[k + 1] - 1;
            if (ist <= ien) {
                if (ann_introns.find(AnnIntKey{csid, ist, ien, st}) == ann_introns.end()) {
                    q_pass_intron[q] = 0;
                    break;
                }
            }
        }
    }

    std::vector<int32_t> comp_cnt(n_q, 0);
    std::vector<uint8_t> eq_flag(n_q, 0);

    int act_thr = std::min(n_threads, 16);
    if (n_q < 512) act_thr = 1;
    std::vector<std::thread> workers;
    workers.reserve(act_thr);
    for (int t = 0; t < act_thr; ++t) {
        int32_t beg = static_cast<int32_t>((static_cast<int64_t>(n_q) * t) / act_thr);
        int32_t end = static_cast<int32_t>((static_cast<int64_t>(n_q) * (t + 1)) / act_thr);
        workers.emplace_back([&, beg, end]() {
            std::vector<int32_t> cands;
            cands.reserve(64);
            for (int32_t q = beg; q < end; ++q) {
                if (!q_pass_intron[q]) continue;
                int32_t qb = q_ptr[q];
                int32_t qe = q_ptr[q + 1];
                int32_t n_qe = qe - qb;
                int32_t csid = q_csr_sid[q];
                int8_t qst = q_strand[q];
                cands.clear();
                csr->query_buckets(csr->cut_span_index, csid, qst, false,
                                   q_cut_mins[q], q_cut_maxe[q],
                                   [&](int32_t tx) { cands.push_back(tx); });
                int32_t c_cnt = 0;
                uint8_t is_eq = 0;
                for (int32_t tx : cands) {
                    int32_t tb = csr->tx_exon_ptr[tx];
                    int32_t te = csr->tx_exon_ptr[tx + 1];
                    int32_t n_te = te - tb;
                    if (n_qe > n_te) continue; // each query cut exon must lie within a distinct subject exon (intron cannot hit subject exon!)

                    // Find first subject exon covering query exon 0
                    int32_t k_q = qb;
                    int32_t k_t = tb;
                    // Advance k_t to cover q_cut[k_q]
                    bool all_within = true;
                    while (k_q < qe && k_t < te) {
                        if (csr->cut_exon_end[k_t] < q_cut_s[k_q]) {
                            k_t++;
                        } else if (csr->cut_exon_start[k_t] <= q_cut_s[k_q] &&
                                   csr->cut_exon_end[k_t] >= q_cut_e[k_q]) {
                            k_q++;
                            k_t++;
                        } else {
                            all_within = false;
                            break;
                        }
                    }
                    if (!all_within || k_q < qe) continue;

                    // Verify no subject cut exon intersects any query intron (myCompatibleTranscription)
                    bool splc_ok = true;
                    if (n_qe >= 2) {
                        int32_t iq = qb;
                        int32_t it = tb;
                        while (iq + 1 < qe && it < te) {
                            int32_t is = q_cut_e[iq] + 1;
                            int32_t ie = q_cut_s[iq + 1] - 1;
                            if (csr->cut_exon_end[it] < is) {
                                it++;
                            } else if (csr->cut_exon_start[it] > ie) {
                                iq++;
                            } else {
                                splc_ok = false;
                                break;
                            }
                        }
                    }
                    if (!splc_ok) continue;

                    c_cnt++;
                    if (n_qe == n_te) {
                        bool exact = true;
                        for (int32_t i = 0; i < n_qe; ++i) {
                            if (q_cut_s[qb + i] != csr->cut_exon_start[tb + i] ||
                                q_cut_e[qb + i] != csr->cut_exon_end[tb + i]) {
                                exact = false;
                                break;
                            }
                        }
                        if (exact) is_eq = 1;
                    }
                }
                comp_cnt[q] = c_cnt;
                eq_flag[q] = is_eq;
            }
        });
    }
    for (auto &th : workers) th.join();

    for (int32_t q = 0; q < n_q; ++q) {
        out_comp[q] = comp_cnt[q];
        out_equal[q] = (eq_flag[q] != 0);
    }
    return Rcpp::List::create(Rcpp::Named("compatible") = out_comp, Rcpp::Named("equal") = out_equal);
    END_RCPP
}

// Fast C++ countPolyATerminals via htslib/faidx.h
extern "C" SEXP bambu_next_count_polya_faidx_cpp(
    SEXP fa_path_sexp,
    SEXP fai_path_sexp,
    SEXP seq_names_sexp,
    SEXP seq_lengths_sexp,
    SEXP start_seq_1b_sexp,
    SEXP start_s_sexp,
    SEXP start_e_sexp,
    SEXP end_seq_1b_sexp,
    SEXP end_s_sexp,
    SEXP end_e_sexp,
    SEXP strand_code_sexp,
    SEXP n_threads_sexp
) {
    BEGIN_RCPP
    std::string fa_path = Rcpp::as<std::string>(fa_path_sexp);
    std::string fai_path = Rf_isNull(fai_path_sexp) ? "" : Rcpp::as<std::string>(fai_path_sexp);
    Rcpp::CharacterVector snames_r(seq_names_sexp);
    Rcpp::IntegerVector slens_r(seq_lengths_sexp);
    Rcpp::IntegerVector s_seq(start_seq_1b_sexp);
    Rcpp::IntegerVector s_s(start_s_sexp);
    Rcpp::IntegerVector s_e(start_e_sexp);
    Rcpp::IntegerVector e_seq(end_seq_1b_sexp);
    Rcpp::IntegerVector e_s(end_s_sexp);
    Rcpp::IntegerVector e_e(end_e_sexp);
    Rcpp::IntegerVector st_code(strand_code_sexp);
    int n_threads = std::max(1, Rcpp::as<int>(n_threads_sexp));

    int32_t n = static_cast<int32_t>(st_code.size());
    Rcpp::IntegerVector nAs(n), nAe(n), nTs(n), nTe(n);
    if (n == 0) {
        return Rcpp::List::create(
            Rcpp::Named("numAstart") = nAs, Rcpp::Named("numAend") = nAe,
            Rcpp::Named("numTstart") = nTs, Rcpp::Named("numTend") = nTe);
    }

    std::vector<std::string> snames(snames_r.size());
    std::vector<int32_t> slens(slens_r.size());
    for (R_xlen_t i = 0; i < snames_r.size(); ++i) snames[i] = Rcpp::as<std::string>(snames_r[i]);
    for (R_xlen_t i = 0; i < slens_r.size(); ++i) slens[i] = (slens_r[i] == NA_INTEGER) ? -1 : slens_r[i];

    // Copy inputs out of SEXPs for parallel worker safety
    std::vector<int32_t> vs_seq(n), vs_s(n), vs_e(n), ve_seq(n), ve_s(n), ve_e(n);
    std::vector<int8_t> vst(n);
    for (int32_t i = 0; i < n; ++i) {
        vs_seq[i] = s_seq[i] - 1;
        // resize(..., width=10, fix='start'/'end', ignore.strand=FALSE)
        int8_t st = static_cast<int8_t>(st_code[i]);
        vst[i] = st;
        if (st == -1) {
            vs_s[i] = std::max(1, s_e[i] - 9);
            vs_e[i] = s_e[i];
            ve_seq[i] = e_seq[i] - 1;
            ve_s[i] = e_s[i];
            ve_e[i] = e_s[i] + 9;
        } else {
            vs_s[i] = s_s[i];
            vs_e[i] = s_s[i] + 9;
            ve_seq[i] = e_seq[i] - 1;
            ve_s[i] = std::max(1, e_e[i] - 9);
            ve_e[i] = e_e[i];
        }
    }

    std::vector<int32_t> outA_s(n, 0), outA_e(n, 0), outT_s(n, 0), outT_e(n, 0);
    int act_thr = std::min(n_threads, 8);
    if (n < 1024) act_thr = 1;
    std::vector<std::thread> workers;
    workers.reserve(act_thr);
    std::vector<std::string> err(act_thr);

    auto count_at_range = [](faidx_t *fai, const char *chr, int32_t pos1b_beg, int32_t pos1b_end, int32_t clen, int32_t &ca, int32_t &ct) {
        ca = 0; ct = 0;
        if (pos1b_beg < 1) pos1b_beg = 1;
        if (clen > 0 && pos1b_end > clen) pos1b_end = clen;
        if (pos1b_beg > pos1b_end) return;
        hts_pos_t len = 0;
        char *seq = faidx_fetch_seq64(fai, chr, static_cast<hts_pos_t>(pos1b_beg - 1),
                                      static_cast<hts_pos_t>(pos1b_end - 1), &len);
        if (seq) {
            for (hts_pos_t k = 0; k < len; ++k) {
                char c = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[k])));
                if (c == 'A') ca++;
                else if (c == 'T') ct++;
            }
            free(seq);
        }
    };

    for (int t = 0; t < act_thr; ++t) {
        int32_t beg = static_cast<int32_t>((static_cast<int64_t>(n) * t) / act_thr);
        int32_t end = static_cast<int32_t>((static_cast<int64_t>(n) * (t + 1)) / act_thr);
        workers.emplace_back([&, t, beg, end]() {
            if (beg >= end) return;
            faidx_t *fai = fai_path.empty()
                ? fai_load3(fa_path.c_str(), nullptr, nullptr, 0)
                : fai_load3(fa_path.c_str(), fai_path.c_str(), nullptr, 0);
            if (!fai) {
                err[t] = "faidx_load3 failed in countPolyATerminals";
                return;
            }
            for (int32_t i = beg; i < end; ++i) {
                int32_t c1 = vs_seq[i];
                int32_t c2 = ve_seq[i];
                int32_t a1 = 0, t1 = 0, a2 = 0, t2 = 0;
                if (c1 >= 0 && c1 < static_cast<int32_t>(snames.size())) {
                    count_at_range(fai, snames[c1].c_str(), vs_s[i], vs_e[i], slens[c1], a1, t1);
                }
                if (c2 >= 0 && c2 < static_cast<int32_t>(snames.size())) {
                    count_at_range(fai, snames[c2].c_str(), ve_s[i], ve_e[i], slens[c2], a2, t2);
                }
                // Wait: for '-' strand, BSgenome::getSeq reverse-complements (swapping A <-> T)
                if (vst[i] == -1) {
                    outA_s[i] = t1; outT_s[i] = a1;
                    outA_e[i] = t2; outT_e[i] = a2;
                } else {
                    outA_s[i] = a1; outT_s[i] = t1;
                    outA_e[i] = a2; outT_e[i] = t2;
                }
            }
            fai_destroy(fai);
        });
    }
    for (auto &th : workers) th.join();
    for (const auto &e : err) {
        if (!e.empty()) Rcpp::stop("FASTA polyA error: %s", e.c_str());
    }

    std::memcpy(nAs.begin(), outA_s.data(), n * sizeof(int32_t));
    std::memcpy(nAe.begin(), outA_e.data(), n * sizeof(int32_t));
    std::memcpy(nTs.begin(), outT_s.data(), n * sizeof(int32_t));
    std::memcpy(nTe.begin(), outT_e.data(), n * sizeof(int32_t));
    return Rcpp::List::create(
        Rcpp::Named("numAstart") = nAs, Rcpp::Named("numAend") = nAe,
        Rcpp::Named("numTstart") = nTs, Rcpp::Named("numTend") = nTe);
    END_RCPP
}
