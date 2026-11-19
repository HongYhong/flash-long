#include <Rcpp.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "stage1/stage1_types.hpp"

namespace bambu_next {
namespace {

struct JuncKey {
    int32_t seq_id;
    int32_t start;
    int32_t end;
    bool operator==(const JuncKey &o) const {
        return seq_id == o.seq_id && start == o.start && end == o.end;
    }
};

struct JuncKeyHash {
    size_t operator()(const JuncKey &k) const {
        uint64_t x = static_cast<uint64_t>(static_cast<uint32_t>(k.seq_id));
        x = (x * 0x9e3779b97f4a7c15ULL) ^ static_cast<uint32_t>(k.start);
        x = (x * 0xbf58476d1ce4e5b9ULL) ^ static_cast<uint32_t>(k.end);
        x ^= (x >> 30);
        return static_cast<size_t>(x);
    }
};

struct LocalJuncRec {
    int32_t seq_id;
    int32_t start;
    int32_t end;
    int32_t plus_score;
    int32_t minus_score;
};

struct ChunkTask {
    int32_t tid;
    int64_t beg;
    int64_t end; // exclusive on core.pos
};

struct WorkerResult {
    std::vector<UnsplicedRead> unspliced;
    std::vector<int32_t> sp_seq_id;
    std::vector<int32_t> sp_start;
    std::vector<int32_t> sp_end;
    std::vector<int32_t> sp_read_id;
    std::vector<int8_t> sp_pos_strand;
    std::vector<int32_t> sp_junc_ptr;
    std::vector<int32_t> sp_local_junc_ids;

    std::vector<LocalJuncRec> local_juncs;
    std::unordered_map<JuncKey, int32_t, JuncKeyHash> junc_map;

    std::vector<std::string> tracked_names;
    std::vector<int32_t> tracked_ids;

    int32_t assigned_read_count = 0; // reads passing seqlevel filter (including badReads/0-exon before badReads removal)
    int64_t n_bad_reads = 0;
};

static inline void extract_exons_from_bam1(
    const bam1_t *b,
    std::vector<std::pair<int32_t, int32_t>> &exons_out
) {
    exons_out.clear();
    const uint32_t *cigar = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;
    int32_t cur_pos = static_cast<int32_t>(b->core.pos) + 1; // 1-based inclusive
    int32_t cur_start = -1;
    int32_t cur_end = -1;

    for (uint32_t i = 0; i < n_cigar; ++i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);
        if (len == 0) continue;
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (cur_start < 0) cur_start = cur_pos;
            cur_end = cur_pos + static_cast<int32_t>(len) - 1;
            cur_pos += static_cast<int32_t>(len);
        } else if (op == BAM_CREF_SKIP) {
            if (cur_start >= 0 && cur_end > cur_start) { // width > 1 filter matching prepareDataFromBam
                exons_out.emplace_back(cur_start, cur_end);
            }
            cur_pos += static_cast<int32_t>(len);
            cur_start = -1;
            cur_end = -1;
        }
        // BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD do not consume reference
    }
    if (cur_start >= 0 && cur_end > cur_start) {
        exons_out.emplace_back(cur_start, cur_end);
    }
}

static inline void process_alignment_record(
    const bam1_t *b,
    const std::vector<int32_t> &tid_to_seq_id,
    const std::vector<int32_t> &seq_id_to_genome_len,
    bool track_reads,
    std::vector<std::pair<int32_t, int32_t>> &exons_buf,
    WorkerResult &w
) {
    if (b->core.tid < 0 || b->core.tid >= static_cast<int32_t>(tid_to_seq_id.size())) return;
    if ((b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY)) != 0) return;

    int32_t seq_id = tid_to_seq_id[b->core.tid];
    if (seq_id < 0) return;

    // Assign local 0-based sequential read ID matching prepareDataFromBam + keepSeqlevels
    int32_t local_rid = w.assigned_read_count++;

    extract_exons_from_bam1(b, exons_buf);
    if (exons_buf.empty()) {
        // Empty after width>1 removal: counted toward read_id (same as GRangesList slot), skipped from badReads
        if (track_reads) {
            w.tracked_names.emplace_back(bam_get_qname(b));
            w.tracked_ids.push_back(local_rid);
        }
        return;
    }

    int32_t read_start = exons_buf.front().first;
    int32_t read_end = exons_buf.back().second;
    int32_t chr_len = seq_id_to_genome_len[seq_id];
    if (chr_len > 0 && read_end >= chr_len) { // matches badReads >= seqlengths(genomeSequence)[chr]
        w.n_bad_reads++;
        return;
    }

    if (track_reads) {
        w.tracked_names.emplace_back(bam_get_qname(b));
        w.tracked_ids.push_back(local_rid);
    }

    int8_t is_pos = (b->core.flag & BAM_FREVERSE) ? 0 : 1;
    if (exons_buf.size() == 1) {
        w.unspliced.push_back(UnsplicedRead{seq_id, read_start, read_end, local_rid, is_pos});
        return;
    }

    w.sp_seq_id.push_back(seq_id);
    w.sp_start.push_back(read_start);
    w.sp_end.push_back(read_end);
    w.sp_read_id.push_back(local_rid);
    w.sp_pos_strand.push_back(is_pos);
    if (w.sp_junc_ptr.empty()) w.sp_junc_ptr.push_back(0);

    size_t n_ex = exons_buf.size();
    for (size_t i = 0; i + 1 < n_ex; ++i) {
        int32_t j_start = exons_buf[i].second + 1;
        int32_t j_end = exons_buf[i + 1].first - 1;
        if (j_start > j_end) continue;
        JuncKey key{seq_id, j_start, j_end};
        auto it = w.junc_map.find(key);
        int32_t jidx;
        if (it == w.junc_map.end()) {
            jidx = static_cast<int32_t>(w.local_juncs.size());
            w.junc_map.emplace(key, jidx);
            w.local_juncs.push_back(LocalJuncRec{seq_id, j_start, j_end, is_pos ? 1 : 0, is_pos ? 0 : 1});
        } else {
            jidx = it->second;
            if (is_pos) w.local_juncs[jidx].plus_score++;
            else w.local_juncs[jidx].minus_score++;
        }
        w.sp_local_junc_ids.push_back(jidx);
    }
    w.sp_junc_ptr.push_back(static_cast<int32_t>(w.sp_local_junc_ids.size()));
}

static bool is_bam_coord_sorted(sam_hdr_t *hdr) {
    kstring_t ks = {0, 0, nullptr};
    bool sorted = false;
    if (sam_hdr_find_tag_hd(hdr, "SO", &ks) == 0 && ks.s != nullptr) {
        sorted = (std::strcmp(ks.s, "coordinate") == 0);
    }
    if (ks.s) free(ks.s);
    return sorted;
}

} // namespace

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
) {
    out = BamSampleData{};
    out.run_name = run_name;
    out.track_reads = track_reads;
    out.ref_seq_names = genome_seqnames;
    out.genome_seq_lengths = genome_seqlengths;

    samFile *fp = sam_open(bam_path.c_str(), "r");
    if (!fp) Rcpp::stop("Failed to open BAM file: %s", bam_path.c_str());
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) {
        sam_close(fp);
        Rcpp::stop("Failed to read BAM header: %s", bam_path.c_str());
    }

    int32_t n_targets = hdr->n_targets;
    std::vector<std::string> bam_sq(n_targets);
    std::unordered_map<std::string, int32_t> genome_map;
    genome_map.reserve(genome_seqnames.size() * 2 + 16);
    for (size_t i = 0; i < genome_seqnames.size(); ++i) {
        genome_map.emplace(genome_seqnames[i], static_cast<int32_t>(i));
    }
    std::unordered_map<std::string, int32_t> ann_map;
    ann_map.reserve(ann_seqnames.size() * 2 + 16);
    for (size_t i = 0; i < ann_seqnames.size(); ++i) {
        ann_map.emplace(ann_seqnames[i], static_cast<int32_t>(i));
    }
    std::unordered_map<std::string, int32_t> bam_sq_map;
    bam_sq_map.reserve(n_targets * 2 + 16);
    bool all_bam_in_ref = true;
    for (int32_t i = 0; i < n_targets; ++i) {
        bam_sq[i] = hdr->target_name[i] ? std::string(hdr->target_name[i]) : "";
        bam_sq_map.emplace(bam_sq[i], i);
        if (genome_map.find(bam_sq[i]) == genome_map.end()) all_bam_in_ref = false;
    }

    // Replicate seqlevelCheckReadsAnnotation(readGrgList, annotations)
    {
        size_t inter_cnt = 0;
        bool all_bam_in_ann = true;
        for (const auto &s : bam_sq) {
            if (ann_map.find(s) != ann_map.end()) inter_cnt++;
            else all_bam_in_ann = false;
        }
        if (inter_cnt == 0) {
            out.warnings.emplace_back(
                "no annotations with matching seqlevel styles, all missing chromosomes will use de-novo annotations");
        }
        if (!all_bam_in_ann) {
            out.warnings.emplace_back(
                "not all chromosomes present in reference annotations, annotations might be incomplete. Please compare objects on the same reference");
        }
    }

    // Determine active seq_names and BAM tid -> seq_id mapping matching keepSeqlevels
    out.all_bam_sq_in_ref = all_bam_in_ref;
    std::vector<int32_t> tid_to_seq_id(n_targets, -1);
    std::vector<int32_t> seq_id_to_genome_len;
    if (all_bam_in_ref) {
        out.seq_names = bam_sq;
        seq_id_to_genome_len.resize(n_targets, -1);
        for (int32_t i = 0; i < n_targets; ++i) {
            tid_to_seq_id[i] = i;
            auto git = genome_map.find(bam_sq[i]);
            if (git != genome_map.end() && git->second < static_cast<int32_t>(genome_seqlengths.size())) {
                seq_id_to_genome_len[i] = genome_seqlengths[git->second];
            }
        }
    } else {
        // intersect(seqlevels(genomeSequence), seqlevels(readGrgList))
        std::vector<std::string> ref_levels;
        for (const auto &gs : genome_seqnames) {
            if (bam_sq_map.find(gs) != bam_sq_map.end()) ref_levels.push_back(gs);
        }
        if (!ann_seqnames.empty()) {
            std::unordered_map<std::string, int> rl_map;
            for (size_t i = 0; i < ref_levels.size(); ++i) rl_map.emplace(ref_levels[i], 1);
            bool all_ann_in = true;
            for (const auto &as : ann_seqnames) {
                if (rl_map.find(as) == rl_map.end()) {
                    all_ann_in = false;
                    break;
                }
            }
            if (!all_ann_in) {
                std::vector<std::string> next_rl;
                for (const auto &s : ref_levels) {
                    if (ann_map.find(s) != ann_map.end()) next_rl.push_back(s);
                }
                ref_levels.swap(next_rl);
                out.warnings.emplace_back(
                    "not all chromosomes from annotations present in reference genome sequence, annotations without reference genomic sequence are dropped");
            }
        }
        out.warnings.emplace_back(
            "not all chromosomes from reads present in reference genome sequence, reads without reference chromosome sequence are dropped");
        out.seq_names = ref_levels;
        seq_id_to_genome_len.resize(ref_levels.size(), -1);
        for (size_t sid = 0; sid < ref_levels.size(); ++sid) {
            int32_t tid = bam_sq_map[ref_levels[sid]];
            tid_to_seq_id[tid] = static_cast<int32_t>(sid);
            auto git = genome_map.find(ref_levels[sid]);
            if (git != genome_map.end() && git->second < static_cast<int32_t>(genome_seqlengths.size())) {
                seq_id_to_genome_len[sid] = genome_seqlengths[git->second];
            }
        }
        out.ref_seq_names = ref_levels;
    }

    // Attempt loading index
    hts_idx_t *idx = nullptr;
    if (!bai_path.empty()) {
        idx = sam_index_load2(fp, bam_path.c_str(), bai_path.c_str());
    } else {
        idx = sam_index_load(fp, bam_path.c_str());
    }
    bool coord_sorted = is_bam_coord_sorted(hdr);
    out.is_indexed = (idx != nullptr);

    std::vector<WorkerResult> worker_results;

    if (idx != nullptr && coord_sorted && n_threads > 1) {
        // Build ordered chromosome windows weighted by .bai mapped counts
        struct TidStat {
            int32_t tid;
            int64_t len;
            uint64_t reads;
        };
        std::vector<TidStat> active_tids;
        uint64_t total_mapped = 0;
        for (int32_t tid = 0; tid < n_targets; ++tid) {
            if (tid_to_seq_id[tid] < 0) continue;
            uint64_t mapped = 0, unmapped = 0;
            hts_idx_get_stat(idx, tid, &mapped, &unmapped);
            if (mapped == 0) {
                // Keep small check if stats are 0 in custom index; include 1 whole chunk just in case
                mapped = 1;
            }
            active_tids.push_back(TidStat{tid, static_cast<int64_t>(hdr->target_len[tid]), mapped});
            total_mapped += mapped;
        }

        int target_chunks = std::max(n_threads * 4, 8);
        std::vector<ChunkTask> chunks;
        chunks.reserve(target_chunks + active_tids.size());
        for (const auto &ts : active_tids) {
            int n_tid_chunks = 1;
            if (total_mapped > 0 && ts.len > 1000000LL) {
                double frac = static_cast<double>(ts.reads) / static_cast<double>(total_mapped);
                n_tid_chunks = std::max(1, static_cast<int>(std::round(frac * target_chunks)));
                n_tid_chunks = std::min(n_tid_chunks, 64);
            }
            if (n_tid_chunks <= 1 || ts.len <= 0) {
                chunks.push_back(ChunkTask{ts.tid, 0, std::numeric_limits<int64_t>::max()});
            } else {
                int64_t step = std::max<int64_t>(250000LL, (ts.len + n_tid_chunks - 1) / n_tid_chunks);
                for (int64_t pos = 0; pos < ts.len; pos += step) {
                    int64_t end = (pos + step >= ts.len) ? std::numeric_limits<int64_t>::max() : (pos + step);
                    chunks.push_back(ChunkTask{ts.tid, pos, end});
                    if (end == std::numeric_limits<int64_t>::max()) break;
                }
            }
        }

        int n_chunks = static_cast<int>(chunks.size());
        worker_results.resize(n_chunks);
        int actual_threads = std::min(n_threads, std::max(1, n_chunks));

        // Dynamic work-stealing across ordered chunks via atomic cursor
        std::atomic<int> next_chunk{0};
        std::vector<std::string> thread_errs(actual_threads);
        std::vector<std::thread> workers;
        workers.reserve(actual_threads);

        for (int t = 0; t < actual_threads; ++t) {
            workers.emplace_back([&, t]() {
                samFile *tfp = sam_open(bam_path.c_str(), "r");
                if (!tfp) {
                    thread_errs[t] = "sam_open failed";
                    return;
                }
                hts_idx_t *tidx = !bai_path.empty()
                    ? sam_index_load2(tfp, bam_path.c_str(), bai_path.c_str())
                    : sam_index_load(tfp, bam_path.c_str());
                if (!tidx) {
                    sam_close(tfp);
                    thread_errs[t] = "sam_index_load failed";
                    return;
                }
                bam1_t *b = bam_init1();
                std::vector<std::pair<int32_t, int32_t>> exons_buf;
                exons_buf.reserve(64);

                while (true) {
                    int cid = next_chunk.fetch_add(1, std::memory_order_relaxed);
                    if (cid >= n_chunks) break;
                    const ChunkTask &ct = chunks[cid];
                    WorkerResult &wr = worker_results[cid];
                    wr.sp_junc_ptr.push_back(0);

                    int64_t itr_end = (ct.end == std::numeric_limits<int64_t>::max())
                        ? static_cast<int64_t>(hdr->target_len[ct.tid]) + 100000000LL
                        : ct.end;
                    hts_itr_t *itr = sam_itr_queryi(tidx, ct.tid, ct.beg, itr_end);
                    if (!itr) continue;
                    int r;
                    while ((r = sam_itr_next(tfp, itr, b)) >= 0) {
                        int64_t pos = static_cast<int64_t>(b->core.pos);
                        if (pos < ct.beg || pos >= ct.end) continue;
                        process_alignment_record(b, tid_to_seq_id, seq_id_to_genome_len, track_reads, exons_buf, wr);
                    }
                    hts_itr_destroy(itr);
                }

                bam_destroy1(b);
                hts_idx_destroy(tidx);
                sam_close(tfp);
            });
        }
        for (auto &th : workers) th.join();
        for (const auto &err : thread_errs) {
            if (!err.empty()) {
                hts_idx_destroy(idx);
                sam_hdr_destroy(hdr);
                sam_close(fp);
                Rcpp::stop("Worker BAM error: %s", err.c_str());
            }
        }
    } else {
        // Fallback: single-reader exact BAM order with multithreaded BGZF decompression
        if (n_threads > 1 && fp->format.compression == bgzf && fp->fp.bgzf != nullptr) {
            bgzf_mt(fp->fp.bgzf, std::min(n_threads, 8), 256);
        }
        worker_results.resize(1);
        WorkerResult &wr = worker_results[0];
        wr.sp_junc_ptr.push_back(0);
        bam1_t *b = bam_init1();
        std::vector<std::pair<int32_t, int32_t>> exons_buf;
        exons_buf.reserve(64);
        while (sam_read1(fp, hdr, b) >= 0) {
            process_alignment_record(b, tid_to_seq_id, seq_id_to_genome_len, track_reads, exons_buf, wr);
        }
        bam_destroy1(b);
    }

    if (idx) hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    sam_close(fp);

    // 1. Prefix-sum assigned_read_count across chunks for 1-based BAM read_id & gather totals
    int n_chunks = static_cast<int>(worker_results.size());
    std::vector<int32_t> chunk_id_offsets(n_chunks + 1, 0);
    size_t total_unspliced = 0;
    size_t total_spliced = 0;
    size_t total_sp_juncs = 0;
    size_t total_tracked = 0;
    for (int c = 0; c < n_chunks; ++c) {
        const WorkerResult &wr = worker_results[c];
        chunk_id_offsets[c + 1] = chunk_id_offsets[c] + wr.assigned_read_count;
        total_unspliced += wr.unspliced.size();
        total_spliced += wr.sp_read_id.size();
        total_sp_juncs += wr.sp_local_junc_ids.size();
        total_tracked += wr.tracked_ids.size();
        out.n_bad_reads += wr.n_bad_reads;
    }
    out.n_valid_reads = chunk_id_offsets[n_chunks];

    if (out.n_bad_reads > 0) {
        out.warnings.emplace_back(
            std::to_string(out.n_bad_reads) +
            " reads are mapped outside the provided genomic regions. These reads will be dropped. Check you are using the same genome used for the alignment");
    }

    // 2. Merge unique junctions across chunks & sort by (seq_id, start, end)
    std::unordered_map<JuncKey, int32_t, JuncKeyHash> merged_jmap;
    size_t est_j = 0;
    for (int c = 0; c < n_chunks; ++c) est_j += worker_results[c].local_juncs.size();
    merged_jmap.reserve(est_j / 2 + 128);

    std::vector<UniqueJunctionRec> unsorted_juncs;
    unsorted_juncs.reserve(est_j / 2 + 128);
    std::vector<std::vector<int32_t>> local_to_merged(n_chunks);

    for (int c = 0; c < n_chunks; ++c) {
        const auto &lj = worker_results[c].local_juncs;
        auto &map_vec = local_to_merged[c];
        map_vec.resize(lj.size());
        for (size_t i = 0; i < lj.size(); ++i) {
            JuncKey k{lj[i].seq_id, lj[i].start, lj[i].end};
            auto it = merged_jmap.find(k);
            if (it == merged_jmap.end()) {
                int32_t midx = static_cast<int32_t>(unsorted_juncs.size());
                merged_jmap.emplace(k, midx);
                UniqueJunctionRec rec;
                rec.seq_id = lj[i].seq_id;
                rec.start = lj[i].start;
                rec.end = lj[i].end;
                rec.plus_score = lj[i].plus_score;
                rec.minus_score = lj[i].minus_score;
                rec.score = rec.plus_score + rec.minus_score;
                unsorted_juncs.push_back(rec);
                map_vec[i] = midx;
            } else {
                int32_t midx = it->second;
                unsorted_juncs[midx].plus_score += lj[i].plus_score;
                unsorted_juncs[midx].minus_score += lj[i].minus_score;
                unsorted_juncs[midx].score += (lj[i].plus_score + lj[i].minus_score);
                map_vec[i] = midx;
            }
        }
    }

    int32_t n_uniq_j = static_cast<int32_t>(unsorted_juncs.size());
    std::vector<int32_t> jperm(n_uniq_j);
    for (int32_t i = 0; i < n_uniq_j; ++i) jperm[i] = i;
    std::sort(jperm.begin(), jperm.end(), [&](int32_t a, int32_t b) {
        const auto &ja = unsorted_juncs[a];
        const auto &jb = unsorted_juncs[b];
        if (ja.seq_id != jb.seq_id) return ja.seq_id < jb.seq_id;
        if (ja.start != jb.start) return ja.start < jb.start;
        return ja.end < jb.end;
    });
    std::vector<int32_t> merged_to_sorted(n_uniq_j);
    out.unique_junctions.resize(n_uniq_j);
    for (int32_t new_i = 0; new_i < n_uniq_j; ++new_i) {
        int32_t old_i = jperm[new_i];
        merged_to_sorted[old_i] = new_i;
        out.unique_junctions[new_i] = unsorted_juncs[old_i];
    }

    // 3. Compute startScore (linear sweep) & endScore (per-chr hash) on sorted unique_junctions
    {
        int32_t i = 0;
        while (i < n_uniq_j) {
            int32_t s = i;
            int32_t seq = out.unique_junctions[i].seq_id;
            int32_t st = out.unique_junctions[i].start;
            int32_t sum_s = 0;
            while (i < n_uniq_j && out.unique_junctions[i].seq_id == seq && out.unique_junctions[i].start == st) {
                sum_s += out.unique_junctions[i].score;
                i++;
            }
            for (int32_t k = s; k < i; ++k) out.unique_junctions[k].start_score = sum_s;
        }
        int32_t n_seq = static_cast<int32_t>(out.seq_names.size());
        i = 0;
        std::unordered_map<int32_t, int32_t> end_sum;
        while (i < n_uniq_j) {
            int32_t s = i;
            int32_t seq = out.unique_junctions[i].seq_id;
            end_sum.clear();
            while (i < n_uniq_j && out.unique_junctions[i].seq_id == seq) {
                end_sum[out.unique_junctions[i].end] += out.unique_junctions[i].score;
                i++;
            }
            for (int32_t k = s; k < i; ++k) {
                out.unique_junctions[k].end_score = end_sum[out.unique_junctions[k].end];
            }
        }
        (void)n_seq;
    }

    // 4. Assemble contiguous unspliced, spliced, and tracked read vectors with remapped IDs
    out.unspliced.reserve(total_unspliced);
    out.spliced.n_reads = static_cast<int32_t>(total_spliced);
    out.spliced.seq_id.resize(total_spliced);
    out.spliced.start.resize(total_spliced);
    out.spliced.end.resize(total_spliced);
    out.spliced.read_id.resize(total_spliced);
    out.spliced.pos_strand.resize(total_spliced);
    out.spliced.junc_ptr.resize(total_spliced + 1);
    out.spliced.junc_ids.resize(total_sp_juncs);
    if (track_reads) {
        out.tracked_read_names.reserve(total_tracked);
        out.tracked_read_ids.reserve(total_tracked);
    }

    size_t sp_off = 0;
    size_t j_off = 0;
    for (int c = 0; c < n_chunks; ++c) {
        WorkerResult &wr = worker_results[c];
        int32_t rid_base = chunk_id_offsets[c] + 1; // 1-based R read_id!

        for (auto &ur : wr.unspliced) {
            ur.read_id += rid_base;
            out.unspliced.push_back(ur);
        }
        size_t n_sp = wr.sp_read_id.size();
        for (size_t i = 0; i < n_sp; ++i) {
            out.spliced.seq_id[sp_off + i] = wr.sp_seq_id[i];
            out.spliced.start[sp_off + i] = wr.sp_start[i];
            out.spliced.end[sp_off + i] = wr.sp_end[i];
            out.spliced.read_id[sp_off + i] = wr.sp_read_id[i] + rid_base;
            out.spliced.pos_strand[sp_off + i] = wr.sp_pos_strand[i];
            out.spliced.junc_ptr[sp_off + i] = static_cast<int32_t>(j_off + wr.sp_junc_ptr[i]);
        }
        const auto &l2m = local_to_merged[c];
        size_t n_lj = wr.sp_local_junc_ids.size();
        for (size_t i = 0; i < n_lj; ++i) {
            out.spliced.junc_ids[j_off + i] = merged_to_sorted[l2m[wr.sp_local_junc_ids[i]]];
        }
        sp_off += n_sp;
        j_off += n_lj;

        if (track_reads) {
            for (size_t i = 0; i < wr.tracked_ids.size(); ++i) {
                out.tracked_read_names.push_back(std::move(wr.tracked_names[i]));
                out.tracked_read_ids.push_back(wr.tracked_ids[i] + rid_base);
            }
        }
        // Release chunk memory immediately to minimize peak RSS during merge
        wr = WorkerResult{};
    }
    out.spliced.junc_ptr[total_spliced] = static_cast<int32_t>(j_off);
}

} // namespace bambu_next
