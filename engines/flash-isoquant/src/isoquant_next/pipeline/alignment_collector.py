"""Locus-streaming AlignmentCollector using `.iqidx` bisect range index and LazyContigView."""

from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from enum import Enum, unique
import logging
import math
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pysam

from isoquant_next.alignment.alignment_info import AlignmentInfo
from isoquant_next.alignment.bam_merger import BAMStreamMerger
from isoquant_next.alignment.exon_corrector import ExonCorrector
from isoquant_next.alignment.illumina_exon_corrector import IlluminaExonCorrector, VoidExonCorrector
from isoquant_next.alignment.read_groups import DefaultReadGrouper
from isoquant_next.annotation.gene_info import GeneInfo, StrandDetector
from isoquant_next.annotation.iqidx import MmapAnnotationIndex
from isoquant_next.assignment.assigner import LongReadAssigner
from isoquant_next.assignment.profiles import CombinedProfileConstructor
from isoquant_next.assignment.types import (
    IsoformMatch,
    MatchClassification,
    MatchEvent,
    MatchEventSubtype,
    PolyAFixer,
    ReadAssignment,
    ReadAssignmentType,
)
from isoquant_next.config.enums import PolyATrimmed
from isoquant_next.core.stats import EnumStats
from isoquant_next.core.utils import interval_len, junctions_from_blocks, overlaps
from isoquant_next.discovery.splice_corrector import extract_relevant_deletions
from isoquant_next.terminal.cage_finder import CagePeakFinder
from isoquant_next.terminal.polya_finder import PolyAFinder

logger = logging.getLogger("IsoQuant")


@unique
class AlignmentType(Enum):
    primary = 1
    secondary = 2
    supplementary = 20
    unaligned = 255


class AbstractAlignmentStorage:
    COVERAGE_BIN = 256

    def __init__(self) -> None:
        self.coverage_dict: Dict[int, int] = defaultdict(int)
        self.region: Optional[Tuple[int, int]] = None

    def reset(self) -> None:
        self.coverage_dict = defaultdict(int)
        self.region = None

    def add_alignment(self, bam_index: int, alignment: pysam.AlignedSegment) -> None:
        bin_start = alignment.reference_start // AbstractAlignmentStorage.COVERAGE_BIN
        alignment_end = (alignment.reference_end or (alignment.reference_start + 1)) - 1
        bin_end = alignment_end // AbstractAlignmentStorage.COVERAGE_BIN
        for i in range(bin_start, bin_end + 1):
            self.coverage_dict[i] += 1
        if not self.region:
            self.region = (alignment.reference_start, alignment_end)
        else:
            self.region = (min(self.region[0], alignment.reference_start), max(self.region[1], alignment_end))

    def get_alignments(self, region: Optional[Tuple[int, int]] = None):
        raise NotImplementedError()

    def get_read_count(self) -> int:
        raise NotImplementedError()

    def alignment_is_not_adjacent(self, alignment: pysam.AlignedSegment) -> bool:
        alignment_end = (alignment.reference_end or (alignment.reference_start + 1)) - 1
        return self.region is not None and not overlaps(self.region, (alignment.reference_start, alignment_end))


class BAMAlignmentStorage(AbstractAlignmentStorage):
    """Low-memory locus buffer that falls back to BAM region fetch for giant coverage loci."""

    def __init__(self, bam_files: Sequence[str], chr_id: str) -> None:
        super().__init__()
        self.bam_files = list(bam_files)
        self.chr_id = chr_id
        self.counter = 0

    def reset(self) -> None:
        super().reset()
        self.counter = 0

    def add_alignment(self, bam_index: int, alignment: pysam.AlignedSegment) -> None:
        super().add_alignment(bam_index, alignment)
        self.counter += 1

    def get_alignments(self, region: Optional[Tuple[int, int]] = None):
        if not region:
            region = self.region
        if region is None:
            return iter(())
        return BAMStreamMerger(
            self.bam_files,
            chr_id=self.chr_id,
            start=max(0, region[0]),
            end=region[1] + 1,
            include_secondary=True,
            include_supplementary=True,
        ).iter_alignments()

    def get_read_count(self) -> int:
        return self.counter


class LocusBufferedAlignmentStorage(AbstractAlignmentStorage):
    """Stores ONLY the currently active non-adjacent locus alignments (O(Locus) RAM, reset per locus)."""

    def __init__(self, bam_files: Sequence[str], chr_id: str, spill_threshold: int = 100_000) -> None:
        super().__init__()
        self.bam_files = list(bam_files)
        self.chr_id = chr_id
        self.spill_threshold = spill_threshold
        self.alignment_start_index: Dict[int, int] = {}
        self.alignment_end_index: Dict[int, int] = {}
        self.counter = 0
        self.alignment_storage: List[Tuple[int, pysam.AlignedSegment]] = []
        self.spilled = False
        self.index_filled = False

    def reset(self) -> None:
        super().reset()
        self.alignment_start_index = {}
        self.alignment_end_index = {}
        self.counter = 0
        self.alignment_storage.clear()
        self.spilled = False
        self.index_filled = False

    def add_alignment(self, bam_index: int, alignment: pysam.AlignedSegment) -> None:
        super().add_alignment(bam_index, alignment)
        if self.spilled:
            self.counter += 1
            return
        if len(self.alignment_storage) >= self.spill_threshold:
            self.spilled = True
            self.alignment_storage.clear()
            self.alignment_start_index.clear()
            self.alignment_end_index.clear()
            self.counter += 1
            return
        bin_start = alignment.reference_start // self.COVERAGE_BIN
        if bin_start not in self.alignment_start_index:
            self.alignment_start_index[bin_start] = self.counter
        alignment_end = (alignment.reference_end or (alignment.reference_start + 1)) - 1
        bin_end = alignment_end // self.COVERAGE_BIN
        if bin_end not in self.alignment_end_index:
            self.alignment_end_index[bin_end] = self.counter
        self.alignment_storage.append((bam_index, alignment))
        self.counter += 1
        self.index_filled = False

    def fill_index(self) -> None:
        if self.index_filled or not self.region:
            return
        cur_idx = len(self.alignment_storage)
        b0 = self.region[0] // self.COVERAGE_BIN
        b1 = self.region[1] // self.COVERAGE_BIN
        for pos in range(b1 + 1, b0 - 1, -1):
            if pos not in self.alignment_start_index:
                self.alignment_start_index[pos] = cur_idx
            else:
                cur_idx = self.alignment_start_index[pos]
        cur_idx = len(self.alignment_storage)
        for pos in range(b1 + 1, b0 - 1, -1):
            if pos not in self.alignment_end_index or self.alignment_end_index[pos] > cur_idx:
                self.alignment_end_index[pos] = cur_idx
            else:
                cur_idx = self.alignment_end_index[pos]
        self.index_filled = True

    def get_alignments(self, region: Optional[Tuple[int, int]] = None):
        if self.spilled:
            target = region or self.region
            if target is None:
                return iter(())
            return BAMStreamMerger(
                self.bam_files,
                chr_id=self.chr_id,
                start=max(0, target[0]),
                end=target[1] + 1,
                include_secondary=True,
                include_supplementary=True,
            ).iter_alignments()
        if not region or region == self.region:
            return iter(self.alignment_storage)
        self.fill_index()
        start_bin = region[0] // self.COVERAGE_BIN
        start_idx = self.alignment_end_index.get(start_bin, 0)
        end_bin = region[1] // self.COVERAGE_BIN
        end_idx = self.alignment_start_index.get(end_bin, len(self.alignment_storage))
        return (
            (bi, aln)
            for bi, aln in self.alignment_storage[start_idx:end_idx]
            if overlaps(region, (aln.reference_start, (aln.reference_end or (aln.reference_start + 1)) - 1))
        )

    def get_read_count(self) -> int:
        return self.counter


class _IqidxChrGeneLookup:
    """Block-indexed O(log G) interval overlap lookup over one chromosome of MmapAnnotationIndex."""

    _BLK = 64

    def __init__(self, iqidx: MmapAnnotationIndex, chr_id: str) -> None:
        self.iqidx = iqidx
        g0, g1 = iqidx.chrom_gene_slice(chr_id)
        self.g0 = g0
        self.g1 = g1
        if g0 < g1:
            self.starts = [int(iqidx.gene_start[i]) for i in range(g0, g1)]
            self.ends = [int(iqidx.gene_end[i]) for i in range(g0, g1)]
            self.block_max = [
                max(self.ends[b : b + self._BLK]) for b in range(0, len(self.ends), self._BLK)
            ]
        else:
            self.starts = []
            self.ends = []
            self.block_max = []

    def gene_indices_in_region(self, start: int, end: int) -> List[int]:
        if not self.starts or end < start:
            return []
        limit = bisect_right(self.starts, end)
        if limit <= 0:
            return []
        hits: List[int] = []
        n_blk = (limit + self._BLK - 1) // self._BLK
        for blk in range(n_blk):
            if self.block_max[blk] < start:
                continue
            lo = blk * self._BLK
            hi = min(lo + self._BLK, limit)
            for i in range(lo, hi):
                if self.ends[i] >= start:
                    hits.append(self.g0 + i)
        return hits


class AlignmentCollector:
    """Streams BAM alignments locus-by-locus against `.iqidx` + `LazyContigView` to yield GeneInfo & assignments."""

    MIN_REGION_LEN = 32768
    MAX_REGION_LEN = 1048576
    MIN_READS_TO_SPLIT = 1024
    ABS_COV_VALLEY = 1
    REL_COV_VALLEY = 0.01

    SMALL_CHR_IDS = ("MT", "chrM", "chrMT")
    SMALL_CHR_LEN = 500_000
    WARN_COVERAGE = 2_000_000

    def __init__(
        self,
        chr_id: str,
        bam_files: Sequence[str],
        params,
        illumina_bam: Optional[Sequence[str]] = None,
        genedb=None,
        chr_record=None,
        read_groupper=None,
        barcode_dict: Optional[Dict[str, Tuple[str, str]]] = None,
        small_chr_max_coverage: int = 1_000_000,
        usual_gene_max_coverage: int = -1,
        string_pools=None,
    ) -> None:
        self.chr_id = chr_id
        self.bam_files = list(bam_files)
        self.params = params
        self.genedb = genedb
        self.chr_record = chr_record
        self.illumina_bam = list(illumina_bam) if illumina_bam else None
        self.string_pools = string_pools
        self._iqidx_lookup = (
            _IqidxChrGeneLookup(genedb, chr_id) if isinstance(genedb, MmapAnnotationIndex) else None
        )

        self.strand_detector = StrandDetector(self.chr_record)
        self.read_groupper = read_groupper if read_groupper is not None else DefaultReadGrouper()
        self.barcode_dict = barcode_dict or {}
        self.barcode_tag = params.barcode_tag if getattr(params, "barcoded_bam", False) else None
        self.umi_tag = params.umi_tag if getattr(params, "barcoded_bam", False) else None
        self.strip_barcode_suffix = getattr(params, "strip_barcode_suffix", False)
        self.small_chr_max_coverage = small_chr_max_coverage
        self.usual_gene_max_coverage = usual_gene_max_coverage
        self.polya_finder = PolyAFinder(self.params.polya_window, self.params.polya_fraction)
        self.polya_fixer = PolyAFixer(self.params)
        cage_file = getattr(params, "cage", None)
        self.cage_finder = (
            CagePeakFinder(cage_file, shift_size=getattr(params, "cage_shift", 50))
            if cage_file
            else None
        )
        self.alignment_stat_counter = EnumStats()

    def process(self):
        bam_merger = BAMStreamMerger(
            self.bam_files,
            chr_id=self.chr_id,
            include_unmapped=False,
            include_secondary=True,
            include_supplementary=True,
        )
        storage = LocusBufferedAlignmentStorage(
            self.bam_files,
            self.chr_id,
            spill_threshold=500_000 if getattr(self.params, "high_memory", False) else 100_000,
        )

        for bam_index, alignment in bam_merger.iter_alignments():
            if alignment.is_secondary:
                self.alignment_stat_counter.add(AlignmentType.secondary)
            elif alignment.is_supplementary:
                self.alignment_stat_counter.add(AlignmentType.supplementary)
            elif alignment.reference_id != -1:
                self.alignment_stat_counter.add(AlignmentType.primary)

            if storage.alignment_is_not_adjacent(alignment):
                assert storage.region is not None
                preceding_genes = self.get_genes_in_region(storage.region)
                aln_end = (alignment.reference_end or (alignment.reference_start + 1)) - 1
                next_genes = self.get_genes_in_region((alignment.reference_start, aln_end))
                if not (preceding_genes & next_genes):
                    yield from self.forward_alignments(storage)
                    storage.reset()

            storage.add_alignment(bam_index, alignment)

        if storage.region:
            yield from self.forward_alignments(storage)
            storage.reset()

    def forward_alignments(self, alignment_storage):
        current_region = alignment_storage.region
        split_regions = AlignmentCollector.split_coverage_regions(current_region, alignment_storage)
        if len(split_regions) == 1:
            max_coverage = split_regions[0][1]
            yield self.process_alignments_in_region(current_region, alignment_storage.get_alignments(), max_coverage)
        else:
            for new_region, max_coverage in split_regions:
                yield self.process_alignments_in_region(
                    new_region, alignment_storage.get_alignments(new_region), max_coverage
                )

    def process_alignments_in_region(self, current_region: Tuple[int, int], alignments, max_coverage: int):
        skip_read_fraction = 1
        chr_len = len(self.chr_record) if self.chr_record is not None else 0
        if self.chr_id in self.SMALL_CHR_IDS or (0 < chr_len < self.SMALL_CHR_LEN):
            coverage_cutoff = self.small_chr_max_coverage
        else:
            coverage_cutoff = self.usual_gene_max_coverage

        if max_coverage > coverage_cutoff > 0:
            skip_read_fraction = math.ceil(max_coverage / coverage_cutoff)

        gene_info = self.get_gene_info_for_region(current_region)
        if gene_info.empty():
            assignments = list(self.process_intergenic(alignments, gene_info, current_region, skip_read_fraction))
        else:
            assignments = list(self.process_genic(alignments, gene_info, current_region, skip_read_fraction))
        return gene_info, assignments

    def _populate_common_fields(self, bam_index: int, alignment, alignment_info: AlignmentInfo, read_assignment: ReadAssignment, gene_info, region: Tuple[int, int]) -> None:
        read_id = alignment.query_name
        if alignment_info.exons_changed:
            read_assignment.add_match_attribute(MatchEvent(MatchEventSubtype.aligned_polya_tail))

        read_assignment.gene_info = gene_info
        read_assignment.polya_info = alignment_info.polya_info
        read_assignment.cage_found = len(alignment_info.cage_hits) > 0
        read_assignment.genomic_region = region
        read_assignment.exons = alignment_info.read_exons
        read_assignment.boundary_deletions = extract_relevant_deletions(
            alignment.cigartuples, alignment_info.read_exons[0][0], alignment_info.read_exons
        )
        if read_assignment.corrected_exons is None:
            read_assignment.corrected_exons = alignment_info.read_exons
            read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)

        if read_id in self.barcode_dict:
            read_assignment.barcode, read_assignment.umi = self.barcode_dict[read_id]
        elif self.barcode_tag and alignment.has_tag(self.barcode_tag):
            barcode = str(alignment.get_tag(self.barcode_tag))
            if self.strip_barcode_suffix and "-" in barcode:
                barcode = barcode.rsplit("-", 1)[0]
            read_assignment.barcode = barcode
            read_assignment.umi = str(alignment.get_tag(self.umi_tag)) if alignment.has_tag(self.umi_tag) else None

        group_ids = self.read_groupper.get_group_id(alignment, read_assignment, self.bam_files[bam_index])
        read_assignment.read_group = group_ids if isinstance(group_ids, list) else [group_ids]
        read_assignment.mapped_strand = "-" if alignment.is_reverse else "+"
        read_assignment.strand = self.get_assignment_strand(read_assignment)
        read_assignment.chr_id = self.chr_id
        read_assignment.multimapper = alignment.is_secondary
        read_assignment.mapping_quality = alignment.mapping_quality

        self.add_artificial_polya(read_assignment)
        pi = alignment_info.polya_info
        read_assignment.polyA_found = (
            pi.external_polya_pos != -1
            or pi.external_polyt_pos != -1
            or pi.internal_polya_pos != -1
            or pi.internal_polyt_pos != -1
        )
        AlignmentCollector.import_bam_tags(alignment, read_assignment, self.params.bam_tags)

    def process_intergenic(self, alignments, gene_info, region: Tuple[int, int], skip_read_fraction: int = 1):
        corrector = (
            IlluminaExonCorrector(self.chr_id, region[0], region[1], self.illumina_bam)
            if self.illumina_bam is not None
            else VoidExonCorrector()
        )
        try:
            counter = 0
            for bam_index, alignment in alignments:
                if (
                    alignment.reference_id == -1
                    or alignment.is_supplementary
                    or (not self.params.use_secondary and alignment.is_secondary)
                ):
                    continue
                if self.params.min_mapq and alignment.mapping_quality < self.params.min_mapq:
                    continue

                alignment_info = AlignmentInfo(alignment)
                if not alignment_info.read_exons:
                    continue
                if len(alignment_info.read_exons) <= 2 and (
                    alignment.is_secondary or alignment.mapping_quality < self.params.simple_alignments_mapq_cutoff
                ):
                    continue

                counter += 1
                if skip_read_fraction > 1 and counter % skip_read_fraction != 0:
                    continue

                if self.params.polya_trimmed == PolyATrimmed.none:
                    alignment_info.add_polya_info(self.polya_finder, self.polya_fixer)
                if self.cage_finder is not None:
                    alignment_info.add_cage_info(self.cage_finder)

                read_assignment = ReadAssignment(
                    alignment.query_name,
                    ReadAssignmentType.intergenic,
                    self.string_pools,
                    match=IsoformMatch(MatchClassification.intergenic, string_pools=self.string_pools),
                )
                read_assignment.corrected_exons = corrector.correct_read(alignment_info)
                read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)
                self._populate_common_fields(bam_index, alignment, alignment_info, read_assignment, gene_info, region)
                yield read_assignment
        finally:
            if hasattr(corrector, "close"):
                corrector.close()

    def process_genic(self, alignments, gene_info, region: Tuple[int, int], skip_read_fraction: int = 1):
        assigner = LongReadAssigner(gene_info, self.params, self.string_pools)
        profile_constructor = CombinedProfileConstructor(gene_info, self.params)
        exon_corrector = ExonCorrector(gene_info, self.params, self.chr_record)

        counter = 0
        for bam_index, alignment in alignments:
            if (
                alignment.reference_id == -1
                or alignment.is_supplementary
                or (not self.params.use_secondary and alignment.is_secondary)
            ):
                continue
            if self.params.min_mapq and alignment.mapping_quality < self.params.min_mapq:
                continue

            counter += 1
            if skip_read_fraction > 1 and counter % skip_read_fraction != 0:
                continue

            read_id = alignment.query_name
            alignment_info = AlignmentInfo(alignment)
            if not alignment_info.read_exons:
                continue

            alignment_info.add_polya_info(self.polya_finder, self.polya_fixer)
            if self.cage_finder is not None:
                alignment_info.add_cage_info(self.cage_finder)
            alignment_info.construct_profiles(profile_constructor)
            read_assignment = assigner.assign_to_isoform(read_id, alignment_info.combined_profile)

            if (
                read_assignment.assignment_type
                not in (
                    ReadAssignmentType.unique,
                    ReadAssignmentType.unique_minor_difference,
                    ReadAssignmentType.ambiguous,
                )
                and alignment.mapping_quality < self.params.inconsistent_mapq_cutoff
            ):
                continue

            read_assignment.corrected_exons = exon_corrector.correct_assigned_read(alignment_info, read_assignment)
            read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)
            self._populate_common_fields(bam_index, alignment, alignment_info, read_assignment, gene_info, region)
            AlignmentCollector.check_antisense(read_assignment)
            if self.params.count_exons:
                read_assignment.exon_gene_profile = alignment_info.combined_profile.read_exon_profile.gene_profile
                read_assignment.intron_gene_profile = alignment_info.combined_profile.read_intron_profile.gene_profile
            yield read_assignment

    def get_assignment_strand(self, read_assignment: ReadAssignment) -> str:
        if read_assignment.isoform_matches and read_assignment.assignment_type in (
            ReadAssignmentType.unique,
            ReadAssignmentType.unique_minor_difference,
        ):
            return read_assignment.isoform_matches[0].transcript_strand

        has_polya = (
            read_assignment.polya_info.external_polya_pos != -1
            or read_assignment.polya_info.internal_polya_pos != -1
        )
        has_polyt = (
            read_assignment.polya_info.external_polyt_pos != -1
            or read_assignment.polya_info.internal_polyt_pos != -1
        )
        if len(read_assignment.exons) == 1:
            if has_polya and not has_polyt:
                return "+"
            if has_polyt and not has_polya:
                return "-"
            return "."
        return self.strand_detector.get_strand(read_assignment.corrected_introns, has_polya, has_polyt)

    def get_genes_in_region(self, current_region: Tuple[int, int]) -> Set[int | str]:
        if not self.genedb:
            return set()
        if self._iqidx_lookup is not None:
            return set(self._iqidx_lookup.gene_indices_in_region(current_region[0], current_region[1]))
        return {
            g.id
            for g in self.genedb.region(
                seqid=self.chr_id, start=current_region[0], end=current_region[1], featuretype="gene"
            )
        }

    def get_gene_info_for_region(self, current_region: Tuple[int, int]) -> GeneInfo:
        s1 = max(1, int(current_region[0]))
        e1 = max(s1, int(current_region[1]))
        if not self.genedb:
            gi = GeneInfo.from_region(self.chr_id, s1, e1, self.params.delta)
            if self.params.needs_reference and self.chr_record is not None:
                gi.set_reference_sequence(s1, min(len(self.chr_record), e1), self.chr_record)
            return gi

        if self._iqidx_lookup is not None:
            gene_indices = self._iqidx_lookup.gene_indices_in_region(current_region[0], current_region[1])
            if not gene_indices:
                gi = GeneInfo.from_region(self.chr_id, s1, e1, self.params.delta)
                if self.params.needs_reference and self.chr_record is not None:
                    gi.set_reference_sequence(s1, min(len(self.chr_record), e1), self.chr_record)
                return gi
            gi = GeneInfo.from_iqidx_cluster(self.genedb, gene_indices, delta=self.params.delta, prepare_profiles=True)
            if self.params.needs_reference and self.chr_record is not None:
                ref_s = max(1, min(gi.start, s1))
                ref_e = min(len(self.chr_record), max(gi.end, e1))
                gi.set_reference_sequence(ref_s, ref_e, self.chr_record)
            return gi

        gene_list = list(
            self.genedb.region(seqid=self.chr_id, start=current_region[0], end=current_region[1], featuretype="gene")
        )
        if not gene_list:
            gi = GeneInfo.from_region(self.chr_id, s1, e1, self.params.delta)
            if self.params.needs_reference and self.chr_record is not None:
                gi.set_reference_sequence(s1, min(len(self.chr_record), e1), self.chr_record)
            return gi
        gene_list.sort(key=lambda x: x.start)
        gi = GeneInfo(gene_list, self.genedb, self.params.delta)
        if self.params.needs_reference and self.chr_record is not None:
            gi.set_reference_sequence(max(1, current_region[0]), min(len(self.chr_record), current_region[1]), self.chr_record)
        return gi

    def add_artificial_polya(self, read_assignment: ReadAssignment) -> None:
        if self.params.polya_trimmed == PolyATrimmed.stranded:
            if read_assignment.strand == "+":
                read_assignment.polya_info.external_polya_pos = read_assignment.corrected_exons[-1][1] + 1
            elif read_assignment.strand == "-":
                read_assignment.polya_info.external_polyt_pos = read_assignment.corrected_exons[0][0] - 1
        elif self.params.polya_trimmed == PolyATrimmed.all:
            if read_assignment.mapped_strand == "+":
                read_assignment.polya_info.external_polya_pos = read_assignment.corrected_exons[-1][1] + 1
            elif read_assignment.mapped_strand == "-":
                read_assignment.polya_info.external_polyt_pos = read_assignment.corrected_exons[0][0] - 1

    @staticmethod
    def split_coverage_regions(genomic_region: Tuple[int, int], alignment_storage) -> List[Tuple[Tuple[int, int], int]]:
        if (
            interval_len(genomic_region) < AlignmentCollector.MIN_REGION_LEN
            or alignment_storage.get_read_count() < AlignmentCollector.MIN_READS_TO_SPLIT
        ):
            max_coverage = max(alignment_storage.coverage_dict.values()) if alignment_storage.coverage_dict else 0
            return [(genomic_region, max_coverage)]

        split_regions = []
        coverage_dict = alignment_storage.coverage_dict
        coverage_positions = sorted(coverage_dict.keys())
        current_start = coverage_positions[0]
        min_bins = int(AlignmentCollector.MIN_REGION_LEN / AbstractAlignmentStorage.COVERAGE_BIN)
        max_bins = int(AlignmentCollector.MAX_REGION_LEN / AbstractAlignmentStorage.COVERAGE_BIN)
        pos = current_start + 1
        max_cov = coverage_dict[current_start]
        abs_cov_valley = AlignmentCollector.ABS_COV_VALLEY
        rel_cov_valley = AlignmentCollector.REL_COV_VALLEY
        prev_upd_pos = pos

        while pos <= coverage_positions[-1]:
            while pos <= coverage_positions[-1] and (
                pos - current_start < min_bins or coverage_dict[pos] > max(abs_cov_valley, max_cov * rel_cov_valley)
            ):
                if pos - prev_upd_pos > max_bins:
                    prev_upd_pos = pos
                    abs_cov_valley += AlignmentCollector.ABS_COV_VALLEY
                    rel_cov_valley += AlignmentCollector.REL_COV_VALLEY
                max_cov = max(max_cov, coverage_dict[pos])
                pos += 1

            split_regions.append(
                (
                    (
                        max(current_start * AbstractAlignmentStorage.COVERAGE_BIN + 1, genomic_region[0]),
                        min(pos * AbstractAlignmentStorage.COVERAGE_BIN, genomic_region[1]),
                    ),
                    max_cov,
                )
            )
            current_start = pos
            prev_upd_pos = pos
            abs_cov_valley = AlignmentCollector.ABS_COV_VALLEY
            rel_cov_valley = AlignmentCollector.REL_COV_VALLEY
            max_cov = coverage_dict[current_start]
            pos = min(current_start + 1, coverage_positions[-1] + 1)

        return split_regions

    @staticmethod
    def check_antisense(read_assignment: ReadAssignment) -> None:
        for match in read_assignment.isoform_matches:
            if match.assigned_transcript is None:
                continue
            if read_assignment.strand == match.transcript_strand or read_assignment.strand == ".":
                continue
            if match.match_classification not in (
                MatchClassification.novel_in_catalog,
                MatchClassification.novel_not_in_catalog,
            ):
                match.match_classification = MatchClassification.antisense
            match.match_subclassifications.append(MatchEvent(MatchEventSubtype.antisense))

    @staticmethod
    def import_bam_tags(bam_record, read_assignment: ReadAssignment, tag_list: Sequence[str]) -> None:
        for tag in tag_list:
            try:
                read_assignment.set_additional_attribute(tag, str(bam_record.get_tag(tag)))
            except KeyError:
                pass
