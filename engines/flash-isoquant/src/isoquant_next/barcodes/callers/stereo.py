"""Stereo-seq and shared-memory spatial barcode detectors (with context cleanup)."""

from __future__ import annotations

import logging
from typing import List, Optional

from isoquant_next.barcodes.common import (
    batch_str_to_2bit_chunked,
    detect_exact_positions,
    detect_first_exact_positions,
    find_candidate_with_max_score_ssw,
    find_polyt_start,
    reverese_complement,
)
from isoquant_next.barcodes.indexers import (
    ArrayKmerIndexer,
    Dict2BitKmerIndexer,
    SharedMemoryArray2BitKmerIndexer,
)
from isoquant_next.barcodes.results import (
    LinkerBarcodeDetectionResult,
    SplittingBarcodeDetectionResult,
    TSOBarcodeDetectionResult,
)

logger = logging.getLogger("IsoQuant")


class StereoBarcodeDetector:
    LINKER = "TTGTCTTCCTAAGAC"
    TSO_PRIMER = "ACTGAGAGGCATGGCGACCTTATCAG"
    PC1_PRIMER = "CTTCCGATCTATGGCGACCTTATCAG"
    BC_LENGTH = 25
    UMI_LEN = 10
    NON_T_UMI_BASES = 0
    UMI_LEN_DELTA = 4
    TERMINAL_MATCH_DELTA = 3
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcodes: List[str]) -> None:
        self.main_primer = self.PC1_PRIMER
        self.pcr_primer_indexer = ArrayKmerIndexer([self.main_primer], kmer_size=6)
        self.linker_indexer = ArrayKmerIndexer([self.LINKER], kmer_size=5)
        self.strict_linker_indexer = ArrayKmerIndexer([self.LINKER], kmer_size=6)
        self.barcode_indexer = None
        if barcodes:
            bit_barcodes = batch_str_to_2bit_chunked(iter(barcodes), seq_len=self.BC_LENGTH)
            self.barcode_indexer = Dict2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        self.min_score = 22

    def close(self) -> None:
        if hasattr(self.barcode_indexer, "close"):
            self.barcode_indexer.close()

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
                l_occ, min_score=12, start_delta=self.TERMINAL_MATCH_DELTA, end_delta=self.TERMINAL_MATCH_DELTA,
            )
        if linker_start is None:
            l_occ = self.strict_linker_indexer.get_occurrences(sequence)
            linker_start, linker_end = detect_exact_positions(
                sequence, 0, len(sequence), self.linker_indexer.k, self.LINKER,
                l_occ, min_score=12, start_delta=self.STRICT_TERMINAL_MATCH_DELTA, end_delta=self.STRICT_TERMINAL_MATCH_DELTA,
            )
        if linker_start is None or linker_end is None:
            return LinkerBarcodeDetectionResult(read_id, polyT=polyt_start)

        if polyt_start == -1:
            presumable = linker_end + self.UMI_LEN
            s0, s1 = presumable - 4, min(len(sequence), presumable + 10)
            p_sub = find_polyt_start(sequence[s0:s1], window_size=5, polya_fraction=1.0)
            if p_sub != -1:
                polyt_start = s0 + p_sub

        p_occ = self.pcr_primer_indexer.get_occurrences(sequence[:linker_start])
        p_start, p_end = detect_exact_positions(
            sequence, 0, linker_start, self.pcr_primer_indexer.k, self.main_primer,
            p_occ, min_score=12, end_delta=self.TERMINAL_MATCH_DELTA,
        )
        if p_start is None or p_end is None:
            p_end = linker_start - self.BC_LENGTH - 1
        if p_end < 0:
            return LinkerBarcodeDetectionResult(
                read_id, polyT=polyt_start, primer=-1, linker_start=linker_start, linker_end=linker_end
            )

        bc_s = p_end + 1
        bc_e = linker_start - 1
        if abs((bc_e - bc_s) - self.BC_LENGTH) > 10:
            return LinkerBarcodeDetectionResult(
                read_id, polyT=polyt_start, primer=p_end, linker_start=linker_start, linker_end=linker_end
            )
        potential_barcode = sequence[bc_s : bc_e + 1]
        assert self.barcode_indexer is not None
        matching = self.barcode_indexer.get_occurrences(potential_barcode, max_hits=10, min_kmers=2)
        barcode, bc_score, _, _ = find_candidate_with_max_score_ssw(
            matching, potential_barcode, min_score=self.min_score, sufficient_score=self.BC_LENGTH - 1
        )
        if barcode is None:
            return LinkerBarcodeDetectionResult(
                read_id, polyT=polyt_start, primer=p_end, linker_start=linker_start, linker_end=linker_end
            )

        u_start = linker_end + 1
        u_end = polyt_start - 1
        if u_end - u_start <= 5:
            u_end = u_start + self.UMI_LEN - 1
        umi = None
        good_umi = False
        if u_start + 2 * self.UMI_LEN > u_end > u_start:
            umi = sequence[u_start : u_end + 1]
            good_umi = abs(len(umi) - self.UMI_LEN) <= self.UMI_LEN_DELTA
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


class SharedMemoryStereoBarcodeDetector(StereoBarcodeDetector):
    MIN_BARCODES_FOR_SHARED_MEM = 1_000_000

    def __init__(self, barcodes: List[str]) -> None:
        super().__init__([])
        bit_barcodes = batch_str_to_2bit_chunked(iter(barcodes), seq_len=self.BC_LENGTH)
        self.barcode_count = len(bit_barcodes)
        self.bit_barcodes = None
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            self.bit_barcodes = bit_barcodes
            self.barcode_indexer = Dict2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        else:
            self.barcode_indexer = SharedMemoryArray2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)

    def __getstate__(self):
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            return (self.min_score, self.barcode_count, self.bit_barcodes)
        return (self.min_score, self.barcode_count, self.barcode_indexer.get_sharable_info())

    def __setstate__(self, state):
        self.min_score = state[0]
        super().__init__([])
        self.bit_barcodes = None
        self.barcode_count = state[1]
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            self.bit_barcodes = state[2]
            self.barcode_indexer = Dict2BitKmerIndexer(self.bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        else:
            self.barcode_indexer = SharedMemoryArray2BitKmerIndexer.from_sharable_info(state[2])


class StereoSplittingBarcodeDetector:
    TSO5 = "CCCGCCTCTCAGTACGTCAGCAG"
    LINKER = "TTGTCTTCCTAAGAC"
    TSO_PRIMER = "ACTGAGAGGCATGGCGACCTTATCAG"
    PC1_PRIMER = "CTTCCGATCTATGGCGACCTTATCAG"
    BC_LENGTH = 25
    UMI_LEN = 10
    NON_T_UMI_BASES = 0
    UMI_LEN_DELTA = 3
    TERMINAL_MATCH_DELTA = 3
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcodes: List[str]) -> None:
        self.main_primer = self.PC1_PRIMER
        self.tso5_indexer = ArrayKmerIndexer([self.TSO5], kmer_size=8)
        self.pcr_primer_indexer = ArrayKmerIndexer([self.main_primer], kmer_size=6)
        self.linker_indexer = ArrayKmerIndexer([self.LINKER], kmer_size=5)
        self.strict_linker_indexer = ArrayKmerIndexer([StereoBarcodeDetector.LINKER], kmer_size=7)
        self.barcode_indexer = None
        if barcodes:
            bit_barcodes = batch_str_to_2bit_chunked(iter(barcodes), seq_len=self.BC_LENGTH)
            self.barcode_indexer = Dict2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        self.min_score = 22

    def close(self) -> None:
        if hasattr(self.barcode_indexer, "close"):
            self.barcode_indexer.close()

    def _scan_dir(self, read_id: str, seq: str, strand: str, out: SplittingBarcodeDetectionResult) -> Optional[TSOBarcodeDetectionResult]:
        r = self._find_barcode_umi_fwd(read_id, seq)
        last_r = r
        prev_start = 0
        while r.polyT != -1:
            r.set_strand(strand)
            out.append(r)
            current_start = (r.tso5 + 15) if r.tso5 != -1 else (r.polyT + 100)
            current_start = max(prev_start + 150, current_start)
            prev_start = current_start
            if len(seq) - current_start < 50:
                break
            r = self._find_barcode_umi_fwd(read_id, seq[current_start:])
            r.update_coordinates(current_start)
            last_r = r
        return last_r

    def find_barcode_umi(self, read_id: str, sequence: str) -> SplittingBarcodeDetectionResult:
        out = SplittingBarcodeDetectionResult(read_id)
        last_fwd = self._scan_dir(read_id, sequence, "+", out)
        self._scan_dir(read_id, reverese_complement(sequence), "-", out)
        if out.empty() and last_fwd is not None:
            out.append(last_fwd)
        out.filter()
        return out

    def find_barcode_umi_single(self, read_id: str, sequence: str) -> TSOBarcodeDetectionResult:
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

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> TSOBarcodeDetectionResult:
        polyt_start = find_polyt_start(sequence)
        linker_start: Optional[int] = None
        linker_end: Optional[int] = None
        tso5_start: Optional[int] = None

        if polyt_start != -1:
            l_occ = self.linker_indexer.get_occurrences(sequence[0 : polyt_start + 1])
            linker_start, linker_end = detect_exact_positions(
                sequence, 0, polyt_start + 1, self.linker_indexer.k, self.LINKER,
                l_occ, min_score=10, start_delta=self.TERMINAL_MATCH_DELTA, end_delta=self.TERMINAL_MATCH_DELTA,
            )
            tso_occ = self.tso5_indexer.get_occurrences(sequence[polyt_start + 1 :])
            tso5_start, _ = detect_first_exact_positions(
                sequence, polyt_start + 1, len(sequence), self.tso5_indexer.k, self.TSO5,
                tso_occ, min_score=15, start_delta=self.TERMINAL_MATCH_DELTA, end_delta=self.TERMINAL_MATCH_DELTA,
            )

        if linker_start is None:
            l_occ = self.strict_linker_indexer.get_occurrences(sequence)
            linker_start, linker_end = detect_first_exact_positions(
                sequence, 0, len(sequence), self.linker_indexer.k, StereoBarcodeDetector.LINKER,
                l_occ, min_score=12, start_delta=self.STRICT_TERMINAL_MATCH_DELTA, end_delta=self.STRICT_TERMINAL_MATCH_DELTA,
            )
        if linker_start is None or linker_end is None:
            return TSOBarcodeDetectionResult(read_id, polyT=polyt_start)

        if polyt_start == -1 or polyt_start < linker_start:
            presumable = linker_end + self.UMI_LEN
            s0, s1 = presumable - 4, min(len(sequence), presumable + 10)
            p_sub = find_polyt_start(sequence[s0:s1], window_size=5, polya_fraction=1.0)
            if p_sub != -1:
                polyt_start = s0 + p_sub
            tso_occ = self.tso5_indexer.get_occurrences(sequence[polyt_start + 1 :])
            tso5_start, _ = detect_first_exact_positions(
                sequence, polyt_start + 1, len(sequence), self.tso5_indexer.k, self.TSO5,
                tso_occ, min_score=15, start_delta=self.TERMINAL_MATCH_DELTA, end_delta=self.TERMINAL_MATCH_DELTA,
            )

        if tso5_start is not None:
            l_occ2 = self.strict_linker_indexer.get_occurrences(sequence[polyt_start + 1 : tso5_start])
            nl_s, _ = detect_exact_positions(
                sequence, polyt_start + 1, tso5_start, self.linker_indexer.k, self.LINKER,
                l_occ2, min_score=12, start_delta=self.STRICT_TERMINAL_MATCH_DELTA, end_delta=self.STRICT_TERMINAL_MATCH_DELTA,
            )
            if nl_s is not None and nl_s != -1 and nl_s - polyt_start > 100:
                tso5_start = nl_s - self.BC_LENGTH - len(self.main_primer) - len(self.TSO5)
        else:
            tso5_start = -1

        p_occ = self.pcr_primer_indexer.get_occurrences(sequence[:linker_start])
        p_start, p_end = detect_exact_positions(
            sequence, 0, linker_start, self.pcr_primer_indexer.k, self.main_primer,
            p_occ, min_score=12, end_delta=self.TERMINAL_MATCH_DELTA,
        )
        if p_start is None or p_end is None:
            p_end = linker_start - self.BC_LENGTH - 1
        if p_end < 0:
            return TSOBarcodeDetectionResult(
                read_id, polyT=polyt_start, primer=-1, linker_start=linker_start, linker_end=linker_end, tso=tso5_start
            )

        bc_s = p_end + 1
        bc_e = linker_start - 1
        if abs((bc_e - bc_s) - self.BC_LENGTH) > 10:
            return TSOBarcodeDetectionResult(
                read_id, polyT=polyt_start, primer=p_end, linker_start=linker_start, linker_end=linker_end, tso=tso5_start
            )

        potential_barcode = sequence[bc_s : bc_e + 1]
        if not self.barcode_indexer:
            return TSOBarcodeDetectionResult(
                read_id, potential_barcode, BC_score=0, polyT=polyt_start, primer=p_end,
                linker_start=linker_start, linker_end=linker_end, tso=tso5_start,
            )

        matching = self.barcode_indexer.get_occurrences(potential_barcode, max_hits=10, min_kmers=2)
        barcode, bc_score, _, _ = find_candidate_with_max_score_ssw(
            matching, potential_barcode, min_score=self.min_score, sufficient_score=self.BC_LENGTH - 1
        )
        if barcode is None:
            return TSOBarcodeDetectionResult(
                read_id, polyT=polyt_start, primer=p_end, linker_start=linker_start, linker_end=linker_end, tso=tso5_start
            )

        u_start = linker_end + 1
        u_end = polyt_start - 1
        if u_end - u_start <= 5:
            u_end = u_start + self.UMI_LEN - 1
        umi = None
        good_umi = False
        if u_start + 2 * self.UMI_LEN > u_end > u_start:
            umi = sequence[u_start : u_end + 1]
            good_umi = abs(len(umi) - self.UMI_LEN) <= self.UMI_LEN_DELTA
        if not umi:
            return TSOBarcodeDetectionResult(
                read_id, barcode, BC_score=bc_score, polyT=polyt_start, primer=p_end,
                linker_start=linker_start, linker_end=linker_end, tso=tso5_start,
            )
        return TSOBarcodeDetectionResult(
            read_id, barcode, umi, bc_score, good_umi, polyT=polyt_start, primer=p_end,
            linker_start=linker_start, linker_end=linker_end, tso=tso5_start,
        )

    @staticmethod
    def result_type():
        return SplittingBarcodeDetectionResult

    @classmethod
    def header(cls) -> str:
        return cls.result_type().header()


class SharedMemoryStereoSplittingBarcodeDetector(StereoSplittingBarcodeDetector):
    MIN_BARCODES_FOR_SHARED_MEM = 1_000_000

    def __init__(self, barcodes: List[str]) -> None:
        super().__init__([])
        bit_barcodes = batch_str_to_2bit_chunked(iter(barcodes), seq_len=self.BC_LENGTH)
        self.barcode_count = len(bit_barcodes)
        self.bit_barcodes = None
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            self.bit_barcodes = bit_barcodes
            self.barcode_indexer = Dict2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        else:
            self.barcode_indexer = SharedMemoryArray2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)

    def __getstate__(self):
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            return (self.min_score, self.barcode_count, self.bit_barcodes)
        return (self.min_score, self.barcode_count, self.barcode_indexer.get_sharable_info())

    def __setstate__(self, state):
        self.min_score = state[0]
        super().__init__([])
        self.bit_barcodes = None
        self.barcode_count = state[1]
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            self.bit_barcodes = state[2]
            self.barcode_indexer = Dict2BitKmerIndexer(self.bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        else:
            self.barcode_indexer = SharedMemoryArray2BitKmerIndexer.from_sharable_info(state[2])
