"""Curio spatial barcode detectors using Cython Smith-Waterman backend."""

from __future__ import annotations

from typing import List, Optional, Tuple, Union

from isoquant_next.barcodes.common import (
    detect_exact_positions,
    find_candidate_with_max_score_ssw,
    find_polyt_start,
    reverese_complement,
    str_to_2bit,
)
from isoquant_next.barcodes.indexers import Array2BitKmerIndexer, ArrayKmerIndexer, KmerIndexer
from isoquant_next.barcodes.results import LinkerBarcodeDetectionResult


class CurioBarcodeDetector:
    LINKER = "TCTTCAGCGTTCCCGAGA"
    PCR_PRIMER = "TACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 8
    RIGHT_BC_LENGTH = 6
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
    UMI_LEN = 9
    NON_T_UMI_BASES = 2
    UMI_LEN_DELTA = 2
    TERMINAL_MATCH_DELTA = 2
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcode_list: Union[List[str], Tuple[List[str], List[str]]]) -> None:
        self.pcr_primer_indexer = ArrayKmerIndexer([self.PCR_PRIMER], kmer_size=6)
        self.linker_indexer = ArrayKmerIndexer([self.LINKER], kmer_size=5)
        if isinstance(barcode_list, tuple):
            joint = [b1 + b2 for b1 in barcode_list[0] for b2 in barcode_list[1]]
        else:
            joint = barcode_list
        if len(joint) <= 100_000:
            self.barcode_indexer = ArrayKmerIndexer(joint, kmer_size=6)
            self.max_barcodes_hits = 20
            self.min_matching_kmers = 1
        else:
            self.barcode_indexer = Array2BitKmerIndexer([str_to_2bit(b) for b in joint], kmer_size=6, seq_len=self.BC_LENGTH)
            self.max_barcodes_hits = 10
            self.min_matching_kmers = 2
        self.min_score = 13

    def find_barcode_umi(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
        r = self._find_barcode_umi_fwd(read_id, sequence)
        if r.polyT != -1:
            r.set_strand("+")
        if r.is_valid():
            return r
        rev = self._find_barcode_umi_fwd(read_id, reverese_complement(sequence))
        if rev.polyT != -1:
            rev.set_strand("-")
        if rev.is_valid():
            return rev
        return r if r.more_informative_than(rev) else rev

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
        polyt_start = find_polyt_start(sequence)
        linker_start: Optional[int] = None
        linker_end: Optional[int] = None
        if polyt_start != -1:
            l_occ = self.linker_indexer.get_occurrences(sequence[0 : polyt_start + 1])
            linker_start, linker_end = detect_exact_positions(
                sequence, 0, polyt_start + 1, self.linker_indexer.k, self.LINKER,
                l_occ, min_score=11, start_delta=self.TERMINAL_MATCH_DELTA, end_delta=self.TERMINAL_MATCH_DELTA,
            )
        if linker_start is None:
            l_occ = self.linker_indexer.get_occurrences(sequence)
            linker_start, linker_end = detect_exact_positions(
                sequence, 0, len(sequence), self.linker_indexer.k, self.LINKER,
                l_occ, min_score=14, start_delta=self.STRICT_TERMINAL_MATCH_DELTA, end_delta=self.STRICT_TERMINAL_MATCH_DELTA,
            )
        if linker_start is None or linker_end is None:
            return LinkerBarcodeDetectionResult(read_id, polyT=polyt_start)

        if polyt_start == -1 or polyt_start - linker_end > self.RIGHT_BC_LENGTH + self.UMI_LEN + 10:
            presumable = linker_end + self.RIGHT_BC_LENGTH + self.UMI_LEN
            s0, s1 = presumable - 4, min(len(sequence), presumable + 10)
            p_sub = find_polyt_start(sequence[s0:s1], window_size=5, polya_fraction=1.0)
            if p_sub != -1:
                polyt_start = s0 + p_sub

        p_occ = self.pcr_primer_indexer.get_occurrences(sequence[:linker_start])
        p_start, p_end = detect_exact_positions(
            sequence, 0, linker_start, self.pcr_primer_indexer.k, self.PCR_PRIMER,
            p_occ, min_score=5, end_delta=self.TERMINAL_MATCH_DELTA,
        )
        if p_start is None or p_end is None:
            p_start, p_end = -1, -1

        bc_start_pos = p_end + 1 if p_start != -1 else linker_start - self.LEFT_BC_LENGTH
        bc_end_pos = linker_end + self.RIGHT_BC_LENGTH + 1
        potential_barcode = sequence[bc_start_pos:linker_start] + sequence[linker_end + 1 : bc_end_pos + 1]
        matching = self.barcode_indexer.get_occurrences(
            potential_barcode, max_hits=self.max_barcodes_hits, min_kmers=self.min_matching_kmers
        )
        barcode, bc_score, _, bc_end = find_candidate_with_max_score_ssw(
            matching, potential_barcode, min_score=self.min_score
        )
        if barcode is None:
            return LinkerBarcodeDetectionResult(
                read_id, polyT=polyt_start, primer=p_end, linker_start=linker_start, linker_end=linker_end
            )

        read_bc_end = bc_start_pos + bc_end - 1 + (linker_end - linker_start + 1)
        u_start = read_bc_end + 1
        u_end = polyt_start - 1
        if u_end - u_start <= 5:
            u_end = u_start + self.UMI_LEN - 1
        umi = sequence[u_start : u_end + 1]
        good_umi = (
            self.UMI_LEN - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LEN + self.UMI_LEN_DELTA
            and all(x != "T" for x in umi[-self.NON_T_UMI_BASES :])
        )
        if not umi:
            return LinkerBarcodeDetectionResult(
                read_id, barcode, BC_score=bc_score, polyT=polyt_start, primer=p_end,
                linker_start=linker_start, linker_end=linker_end,
            )
        return LinkerBarcodeDetectionResult(
            read_id, barcode, umi, bc_score, good_umi, polyT=polyt_start, primer=p_end,
            linker_start=linker_start, linker_end=linker_end,
        )

    @staticmethod
    def result_type():
        return LinkerBarcodeDetectionResult

    @classmethod
    def header(cls) -> str:
        return cls.result_type().header()


class CurioIlluminaDetector:
    LINKER = "TCTTCAGCGTTCCCGAGA"
    PCR_PRIMER = "TACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 8
    RIGHT_BC_LENGTH = 6
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
    MIN_BC_LEN = BC_LENGTH - 4
    UMI_LEN = 9
    NON_T_UMI_BASES = 2
    UMI_LEN_DELTA = 2
    SCORE_DIFF = 1
    TERMINAL_MATCH_DELTA = 1
    STRICT_TERMINAL_MATCH_DELTA = 0

    def __init__(self, joint_barcode_list: List[str]) -> None:
        self.pcr_primer_indexer = KmerIndexer([CurioBarcodeDetector.PCR_PRIMER], kmer_size=6)
        self.linker_indexer = KmerIndexer([CurioBarcodeDetector.LINKER], kmer_size=5)
        if len(joint_barcode_list) <= 100_000:
            self.barcode_indexer = ArrayKmerIndexer(joint_barcode_list, kmer_size=5)
            self.max_barcodes_hits = 20
        else:
            self.barcode_indexer = Array2BitKmerIndexer([str_to_2bit(b) for b in joint_barcode_list], kmer_size=5, seq_len=self.BC_LENGTH)
            self.max_barcodes_hits = 10
        self.min_matching_kmers = 2
        self.min_score = 14

    def find_barcode_umi(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
        r = self._find_barcode_umi_fwd(read_id, sequence)
        if r.polyT != -1:
            r.set_strand("+")
        if r.is_valid():
            return r
        rev = self._find_barcode_umi_fwd(read_id, reverese_complement(sequence))
        if rev.polyT != -1:
            rev.set_strand("-")
        if rev.is_valid():
            return rev
        return r if r.more_informative_than(rev) else rev

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
        l_occ = self.linker_indexer.get_occurrences(sequence)
        linker_start, linker_end = detect_exact_positions(
            sequence, 0, len(sequence), self.linker_indexer.k, self.LINKER,
            l_occ, min_score=14, start_delta=self.TERMINAL_MATCH_DELTA, end_delta=self.TERMINAL_MATCH_DELTA,
        )
        if linker_start is None or linker_end is None:
            return LinkerBarcodeDetectionResult(read_id)

        presumable = linker_end + self.RIGHT_BC_LENGTH + self.UMI_LEN
        s0, s1 = presumable - 4, min(len(sequence), presumable + 10)
        polyt_start = find_polyt_start(sequence[s0:s1], window_size=5, polya_fraction=1.0)
        if polyt_start != -1:
            polyt_start += s0

        bc_s = max(0, linker_start - self.LEFT_BC_LENGTH)
        bc_e = min(len(sequence), linker_end + self.RIGHT_BC_LENGTH + 1)
        potential_barcode = sequence[bc_s:linker_start] + sequence[linker_end + 1 : bc_e]
        if len(potential_barcode) < self.MIN_BC_LEN:
            return LinkerBarcodeDetectionResult(read_id, linker_start=linker_start, linker_end=linker_end)

        matching = self.barcode_indexer.get_occurrences(
            potential_barcode, max_hits=self.max_barcodes_hits, min_kmers=self.min_matching_kmers
        )
        barcode, bc_score, _, bc_end = find_candidate_with_max_score_ssw(
            matching, potential_barcode, min_score=len(potential_barcode) - 1, score_diff=self.SCORE_DIFF
        )
        if barcode is None:
            return LinkerBarcodeDetectionResult(
                read_id, polyT=polyt_start, primer=-1, linker_start=linker_start, linker_end=linker_end
            )

        read_bc_end = bc_s + bc_end - 1 + (linker_end - linker_start + 1)
        u_start = read_bc_end + 1
        u_end = polyt_start - 1
        if polyt_start != -1 or u_end - u_start <= 5:
            u_end = u_start + self.UMI_LEN - 1
        umi = sequence[u_start : min(u_end + 1, len(sequence))]
        good_umi = (
            self.UMI_LEN - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LEN + self.UMI_LEN_DELTA
            and all(x != "T" for x in umi[-self.NON_T_UMI_BASES :])
        )
        if not umi:
            return LinkerBarcodeDetectionResult(
                read_id, barcode, BC_score=bc_score, polyT=polyt_start, primer=-1,
                linker_start=linker_start, linker_end=linker_end,
            )
        return LinkerBarcodeDetectionResult(
            read_id, barcode, umi, bc_score, good_umi, polyT=polyt_start, primer=-1,
            linker_start=linker_start, linker_end=linker_end,
        )

    @staticmethod
    def result_type():
        return LinkerBarcodeDetectionResult

    @classmethod
    def header(cls) -> str:
        return cls.result_type().header()
