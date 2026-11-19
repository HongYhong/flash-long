"""Streaming k-way BAM merger and O(Locus) overlap-free genomic window chunker."""

from __future__ import annotations

import heapq
import logging
from contextlib import ExitStack
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import pysam

from isoquant_next.annotation.iqidx import GeneClusterSlice, MmapAnnotationIndex

logger = logging.getLogger("IsoQuant")


def _open_indexed_bam(bam_path: str) -> pysam.AlignmentFile:
    """Open BAM with automatic .bai/.csi index detection and clean errors."""
    csi_path = bam_path + ".csi"
    bai_path = bam_path + ".bai"
    if not os_path_exists(bai_path) and os_path_exists(csi_path):
        return pysam.AlignmentFile(bam_path, "rb", index_filename=csi_path)
    return pysam.AlignmentFile(bam_path, "rb")


def os_path_exists(path: str) -> bool:
    import os

    return os.path.exists(path)


class BAMStreamMerger:
    """
    Lazy k-way heap merger yielding (file_index, alignment) sorted by genomic start coordinate.
    Holds only O(num_bams) head records in memory instead of loading full chromosomes.
    """

    def __init__(
        self,
        bam_files: Sequence[str],
        chr_id: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        include_unmapped: bool = False,
        include_secondary: bool = False,
        include_supplementary: bool = False,
    ) -> None:
        self.bam_files = list(bam_files)
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.include_unmapped = include_unmapped
        self.include_secondary = include_secondary
        self.include_supplementary = include_supplementary

    def _keep_alignment(self, aln: pysam.AlignedSegment) -> bool:
        if not self.include_unmapped and aln.is_unmapped:
            return False
        if not self.include_secondary and aln.is_secondary:
            return False
        if not self.include_supplementary and aln.is_supplementary:
            return False
        return True

    def iter_alignments(self) -> Iterator[Tuple[int, pysam.AlignedSegment]]:
        with ExitStack() as stack:
            iterators = []
            for file_idx, path in enumerate(self.bam_files):
                af = stack.enter_context(_open_indexed_bam(path))
                if self.chr_id is not None:
                    if self.chr_id not in af.references:
                        continue
                    it = af.fetch(self.chr_id, start=self.start, stop=self.end)
                else:
                    it = af.fetch(until_eof=True)
                iterators.append((file_idx, iter(it)))

            heap: List[Tuple[int, int, int, int, pysam.AlignedSegment, Iterator]] = []
            seq_no = 0
            for file_idx, it in iterators:
                for aln in it:
                    if self._keep_alignment(aln):
                        pos = aln.reference_start if not aln.is_unmapped else (1 << 60)
                        tid = aln.reference_id if not aln.is_unmapped else (1 << 30)
                        heapq.heappush(heap, (tid, pos, seq_no, file_idx, aln, it))
                        seq_no += 1
                        break

            while heap:
                _tid, _pos, _seq, file_idx, aln, it = heapq.heappop(heap)
                yield file_idx, aln
                for nxt in it:
                    if self._keep_alignment(nxt):
                        pos = nxt.reference_start if not nxt.is_unmapped else (1 << 60)
                        tid = nxt.reference_id if not nxt.is_unmapped else (1 << 30)
                        heapq.heappush(heap, (tid, pos, seq_no, file_idx, nxt, it))
                        seq_no += 1
                        break


BAMOnlineMerger = BAMStreamMerger


@dataclass
class LocusWindow:
    """Bounded genomic slice pairing disjoint gene clusters with overlapping BAM alignments."""

    chr_id: str
    start: int
    end: int
    gene_clusters: Tuple[GeneClusterSlice, ...]
    alignments: List[Tuple[int, pysam.AlignedSegment]]
    split_for_density: bool = False


class LocusChunker:
    """
    Stream sorted alignments alongside `.iqidx` gene clusters on a single chromosome,
    splitting into bounded `LocusWindow` batches at clean intergenic/alignment gaps
    while enforcing a high-density watermark (`max_reads_per_locus`) to prevent OOM.
    """

    def __init__(
        self,
        chr_id: str,
        max_reads_per_locus: int = 50_000,
        min_gap: int = 0,
    ) -> None:
        self.chr_id = chr_id
        self.max_reads_per_locus = max(1, int(max_reads_per_locus))
        self.min_gap = max(0, int(min_gap))

    def iter_windows(
        self,
        alignment_stream: Iterable[Tuple[int, pysam.AlignedSegment]],
        iqidx: Optional[MmapAnnotationIndex] = None,
        clusters: Optional[Iterable[GeneClusterSlice]] = None,
    ) -> Iterator[LocusWindow]:
        if clusters is None:
            cluster_iter = iter(iqidx.iter_gene_clusters(self.chr_id)) if iqidx is not None else iter(())
        else:
            cluster_iter = iter(clusters)

        next_cluster: Optional[GeneClusterSlice] = next(cluster_iter, None)
        buf_alns: List[Tuple[int, pysam.AlignedSegment]] = []
        buf_clusters: List[GeneClusterSlice] = []
        win_start = -1
        win_end = -1
        last_split_density = False

        def _flush_standalone_clusters_before(pos: int) -> Iterator[LocusWindow]:
            nonlocal next_cluster
            while next_cluster is not None and next_cluster.end + self.min_gap < pos:
                c = next_cluster
                next_cluster = next(cluster_iter, None)
                yield LocusWindow(self.chr_id, c.start, c.end, (c,), [], split_for_density=False)

        def _absorb_overlapping_clusters(upper_bound: int) -> int:
            nonlocal next_cluster
            cur_end = upper_bound
            while next_cluster is not None and next_cluster.start <= cur_end + self.min_gap:
                c = next_cluster
                next_cluster = next(cluster_iter, None)
                buf_clusters.append(c)
                if c.end > cur_end:
                    cur_end = c.end
            return cur_end

        for file_idx, aln in alignment_stream:
            r_start = int(aln.reference_start) + 1
            r_end = int(aln.reference_end) if aln.reference_end is not None else r_start

            if not buf_alns:
                yield from _flush_standalone_clusters_before(r_start)
                win_start = r_start
                win_end = r_end
                if next_cluster is not None and next_cluster.start <= win_end + self.min_gap:
                    win_start = min(win_start, next_cluster.start)
                win_end = _absorb_overlapping_clusters(win_end)
                buf_alns.append((file_idx, aln))
                continue

            # Check whether this read starts strictly past the active locus window & any overlapping cluster
            has_gap = r_start > win_end + self.min_gap
            density_cap_hit = len(buf_alns) >= self.max_reads_per_locus

            if has_gap or density_cap_hit:
                yield LocusWindow(
                    self.chr_id,
                    win_start,
                    win_end,
                    tuple(buf_clusters),
                    buf_alns,
                    split_for_density=(density_cap_hit and not has_gap),
                )
                last_split_density = density_cap_hit and not has_gap
                buf_alns = []
                # Keep overlapping gene cluster(s) only when forced split occurred mid-cluster
                if last_split_density:
                    buf_clusters = [c for c in buf_clusters if c.end >= r_start]
                    win_start = r_start
                    win_end = _absorb_overlapping_clusters(max([r_end] + [c.end for c in buf_clusters]))
                else:
                    buf_clusters = []
                    yield from _flush_standalone_clusters_before(r_start)
                    win_start = r_start
                    win_end = r_end
                    if next_cluster is not None and next_cluster.start <= win_end + self.min_gap:
                        win_start = min(win_start, next_cluster.start)
                    win_end = _absorb_overlapping_clusters(win_end)
                buf_alns.append((file_idx, aln))
            else:
                buf_alns.append((file_idx, aln))
                if r_end > win_end:
                    win_end = r_end
                win_end = _absorb_overlapping_clusters(win_end)

        if buf_alns:
            yield LocusWindow(
                self.chr_id,
                win_start,
                win_end,
                tuple(buf_clusters),
                buf_alns,
                split_for_density=False,
            )

        while next_cluster is not None:
            c = next_cluster
            next_cluster = next(cluster_iter, None)
            yield LocusWindow(self.chr_id, c.start, c.end, (c,), [], split_for_density=False)
