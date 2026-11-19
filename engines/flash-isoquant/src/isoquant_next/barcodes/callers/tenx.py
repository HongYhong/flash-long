"""10x Genomics (v3, v2, split) and Visium HD spatial barcode detectors."""

from __future__ import annotations

import logging
from typing import List, Optional

from isoquant_next.barcodes.common import (
    detect_exact_positions,
    detect_first_exact_positions,
    find_candidate_with_max_score_ssw,
    find_candidate_with_max_score_ssw_var_len,
    find_optimal_kmer_size,
    find_polyt_start,
    reverese_complement,
    str_to_2bit,
)
from isoquant_next.barcodes.indexers import Array2BitKmerIndexer, ArrayKmerIndexer, KmerIndexer
from isoquant_next.barcodes.results import (
    SplittingBarcodeDetectionResult,
    TenXBarcodeDetectionResult,
    TenXSplitBarcodeDetectionResult,
)

logger = logging.getLogger("IsoQuant")


class TenXBarcodeDetector:
    """10x Genomics v3 single-cell barcode detector."""

    TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT"
    R1 = "CTACACGACGCTCTTCCGATCT"
    BARCODE_LEN_10X = 16
    UMI_LEN = 12
    UMI_LEN_DELTA = 2
    TERMINAL_MATCH_DELTA = 2
    STRICT_TERMINAL_MATCH_DELTA = 1
    BARCODE_SEARCH_MARGIN = 30

    def __init__(self, barcode_list: List[str]) -> None:
        self.r1_indexer = KmerIndexer([self.R1], kmer_size=7)
        bc_length = len(barcode_list[0])
        bc_count = len(barcode_list)
        if bc_count > 1_000_000:
            logger.warning("Large barcode whitelist (%d) may impact speed/RAM", bc_count)

        self.k = find_optimal_kmer_size(bc_length, bc_count)
        if bc_count <= 100_000:
            self.barcode_indexer = ArrayKmerIndexer(barcode_list, kmer_size=self.k)
            self.max_barcodes_hits = 20
            self.min_matching_kmers = 1
            self.min_score = 14
            self.score_diff = 0
        else:
            barcode_bit_list = [str_to_2bit(b) for b in barcode_list]
            self.barcode_indexer = Array2BitKmerIndexer(
                barcode_bit_list, kmer_size=self.k, seq_len=self.BARCODE_LEN_10X
            )
            self.max_barcodes_hits = 10
            self.min_matching_kmers = 2
            self.min_score = 16
            self.score_diff = 1

    def find_barcode_umi(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")
        if read_result.is_valid():
            return read_result

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result

        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    def _find_barcode_near_polyt(
        self, read_id: str, sequence: str, polyt_start: int
    ) -> TenXBarcodeDetectionResult:
        search_start = max(0, polyt_start - self.BARCODE_LEN_10X - self.UMI_LEN - self.BARCODE_SEARCH_MARGIN)
        search_region = sequence[search_start:polyt_start]
        if len(search_region) < self.BARCODE_LEN_10X:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)

        matching = self.barcode_indexer.get_occurrences(
            search_region, max_hits=self.max_barcodes_hits, min_kmers=self.min_matching_kmers
        )
        barcode, bc_score, bc_start, bc_end = find_candidate_with_max_score_ssw(
            matching, search_region, min_score=self.min_score, score_diff=self.score_diff
        )
        if barcode is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)

        abs_bc_end = search_start + bc_end
        potential_umi_start = abs_bc_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN - 1
        potential_umi = sequence[potential_umi_start : potential_umi_end + 1]
        good_umi = self.UMI_LEN - self.UMI_LEN_DELTA <= len(potential_umi) <= self.UMI_LEN + self.UMI_LEN_DELTA
        if not potential_umi:
            return TenXBarcodeDetectionResult(read_id, barcode, BC_score=bc_score, polyT=polyt_start)
        return TenXBarcodeDetectionResult(read_id, barcode, potential_umi, bc_score, good_umi, polyT=polyt_start)

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        polyt_start = find_polyt_start(sequence)
        r1_start: Optional[int] = None
        r1_end: Optional[int] = None

        if polyt_start != -1:
            r1_occ = self.r1_indexer.get_occurrences(sequence[0 : polyt_start + 1])
            r1_start, r1_end = detect_exact_positions(
                sequence, 0, polyt_start + 1, self.r1_indexer.k, self.R1,
                r1_occ, min_score=11, end_delta=self.TERMINAL_MATCH_DELTA,
            )
            if r1_start is None:
                r1_start, r1_end = detect_exact_positions(
                    sequence, 0, polyt_start + 1, self.r1_indexer.k, self.R1,
                    r1_occ, min_score=8, end_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                )

        if r1_start is None:
            r1_occ = self.r1_indexer.get_occurrences(sequence)
            r1_start, r1_end = detect_exact_positions(
                sequence, 0, len(sequence), self.r1_indexer.k, self.R1,
                r1_occ, min_score=18,
                start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                end_delta=self.STRICT_TERMINAL_MATCH_DELTA,
            )

        if r1_start is None:
            if polyt_start != -1:
                return self._find_barcode_near_polyt(read_id, sequence, polyt_start)
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)

        assert r1_end is not None
        if polyt_start == -1 or polyt_start - r1_end > self.BARCODE_LEN_10X + self.UMI_LEN + 10:
            presumable = r1_end + self.BARCODE_LEN_10X + self.UMI_LEN
            s0 = presumable - 4
            s1 = min(len(sequence), presumable + 10)
            p_sub = find_polyt_start(sequence[s0:s1], window_size=5, polya_fraction=1.0)
            if p_sub != -1:
                polyt_start = s0 + p_sub

        barcode_start = r1_end + 1
        barcode_end = r1_end + self.BARCODE_LEN_10X + 1
        potential_barcode = sequence[barcode_start : barcode_end + 1]
        matching = self.barcode_indexer.get_occurrences(
            potential_barcode, max_hits=self.max_barcodes_hits, min_kmers=self.min_matching_kmers
        )
        barcode, bc_score, bc_start, bc_end = find_candidate_with_max_score_ssw(
            matching, potential_barcode, min_score=self.min_score, score_diff=self.score_diff
        )
        if barcode is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start, r1=r1_end)

        read_barcode_end = barcode_start + bc_end - 1
        potential_umi_start = read_barcode_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN - 1
        potential_umi = sequence[potential_umi_start : potential_umi_end + 1]
        good_umi = self.UMI_LEN - self.UMI_LEN_DELTA <= len(potential_umi) <= self.UMI_LEN + self.UMI_LEN_DELTA
        if not potential_umi:
            return TenXBarcodeDetectionResult(read_id, barcode, BC_score=bc_score, polyT=polyt_start, r1=r1_end)
        return TenXBarcodeDetectionResult(read_id, barcode, potential_umi, bc_score, good_umi, polyT=polyt_start, r1=r1_end)

    def find_barcode_umi_no_polya(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        return self.find_barcode_umi(read_id, sequence)

    @staticmethod
    def result_type():
        return TenXBarcodeDetectionResult

    @classmethod
    def header(cls) -> str:
        return cls.result_type().header()


class TenXv2BarcodeDetector(TenXBarcodeDetector):
    UMI_LEN = 10


class VisiumHDBarcodeDetector:
    R1 = "CTACACGACGCTCTTCCGATCT"
    BARCODE1_LEN_VIS = 16
    BARCODE2_LEN_VIS = 15
    SEPARATOR_BASES = 2
    TOTAL_BARCODE_LEN_VIS = BARCODE1_LEN_VIS + BARCODE2_LEN_VIS
    UMI_LEN_VIS = 9
    UMI_LEN_DELTA = 2
    TERMINAL_MATCH_DELTA = 2
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcode_pair_list: List[List[str]]) -> None:
        assert len(barcode_pair_list) == 2
        self.r1_indexer = KmerIndexer([self.R1], kmer_size=7)
        self.part1_list = barcode_pair_list[0]
        self.part2_list = barcode_pair_list[1]
        self.part1_barcode_indexer = ArrayKmerIndexer(self.part1_list, kmer_size=7)
        self.part2_barcode_indexer = ArrayKmerIndexer(self.part2_list, kmer_size=7)
        self.max_barcodes_hits = 20
        self.min_matching_kmers = 2
        self.min_score = 13

    def find_barcode_umi(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")
        if read_result.is_valid():
            return read_result
        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result
        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        polyt_start = find_polyt_start(sequence)
        r1_start: Optional[int] = None
        r1_end: Optional[int] = None
        if polyt_start != -1:
            r1_occ = self.r1_indexer.get_occurrences(sequence[0 : polyt_start + 1])
            r1_start, r1_end = detect_exact_positions(
                sequence, 0, polyt_start + 1, self.r1_indexer.k, self.R1,
                r1_occ, min_score=10, end_delta=self.TERMINAL_MATCH_DELTA,
            )
            if r1_start is None:
                r1_start, r1_end = detect_exact_positions(
                    sequence, 0, polyt_start + 1, self.r1_indexer.k, self.R1,
                    r1_occ, min_score=8, end_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                )
        if r1_start is None:
            r1_occ = self.r1_indexer.get_occurrences(sequence)
            r1_start, r1_end = detect_exact_positions(
                sequence, 0, len(sequence), self.r1_indexer.k, self.R1,
                r1_occ, min_score=17,
                start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                end_delta=self.STRICT_TERMINAL_MATCH_DELTA,
            )
        if r1_end is not None and (
            polyt_start == -1 or polyt_start - r1_end > self.TOTAL_BARCODE_LEN_VIS + self.UMI_LEN_VIS + 10
        ):
            presumable = r1_end + self.TOTAL_BARCODE_LEN_VIS + self.UMI_LEN_VIS
            s0 = presumable - 4
            s1 = min(len(sequence), presumable + 10)
            p_sub = find_polyt_start(sequence[s0:s1], window_size=5, polya_fraction=1.0)
            if p_sub != -1:
                polyt_start = s0 + p_sub

        if polyt_start == -1:
            if r1_start is None or r1_end is None:
                return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)
            u_start = r1_end + 1
            u_end = u_start + self.UMI_LEN_VIS - 1
            potential_umi = sequence[u_start : u_end + 1]

            bc1_s = r1_end + self.UMI_LEN_VIS + 1
            bc1_e = bc1_s + self.BARCODE1_LEN_VIS - 1
            pot1 = sequence[bc1_s : bc1_e + 1]
            m1 = self.part1_barcode_indexer.get_occurrences(pot1, max_hits=self.max_barcodes_hits)
            bc1, s1_sc, _, bc1_end = find_candidate_with_max_score_ssw_var_len(m1, pot1, min_score=self.min_score)

            bc2_s = bc1_s + bc1_end + 1
            bc2_e = bc2_s + self.BARCODE2_LEN_VIS - 1
            pot2 = sequence[bc2_s : bc2_e + 1]
            m2 = self.part2_barcode_indexer.get_occurrences(pot2, max_hits=self.max_barcodes_hits)
            bc2, s2_sc, _, _ = find_candidate_with_max_score_ssw_var_len(m2, pot2, min_score=self.min_score)

            if bc1 is None or bc2 is None:
                return TenXBarcodeDetectionResult(read_id, polyT=polyt_start, r1=r1_end)
            return TenXBarcodeDetectionResult(
                read_id, f"{bc1}|{bc2}", potential_umi, s1_sc + s2_sc,
                UMI_good=True, polyT=polyt_start, r1=r1_end,
            )

        bc2_e = polyt_start - 1 - self.SEPARATOR_BASES
        bc2_s = bc2_e - self.BARCODE2_LEN_VIS + 1
        pot2 = sequence[bc2_s : bc2_e + 1]
        m2 = self.part2_barcode_indexer.get_occurrences(pot2, max_hits=self.max_barcodes_hits)
        bc2, s2_sc, bc2_start, _ = find_candidate_with_max_score_ssw_var_len(m2, pot2, min_score=self.min_score)

        bc1_e = bc2_s + bc2_start - 1
        bc1_s = bc1_e - self.BARCODE1_LEN_VIS + 1
        pot1 = sequence[bc1_s : bc1_e + 1]
        m1 = self.part1_barcode_indexer.get_occurrences(pot1, max_hits=self.max_barcodes_hits)
        bc1, s1_sc, bc1_start, _ = find_candidate_with_max_score_ssw_var_len(m1, pot1, min_score=self.min_score)

        u_end = bc1_s + bc1_start - 1
        u_start = (r1_end + 1) if r1_end is not None else max(0, u_end - self.UMI_LEN_VIS)
        umi_good = abs(u_end - u_start + 1 - self.UMI_LEN_VIS) <= self.UMI_LEN_DELTA
        potential_umi = sequence[u_start : u_end + 1]
        r1_val = r1_end if r1_end is not None else -1
        if bc1 is None or bc2 is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start, r1=r1_val)
        return TenXBarcodeDetectionResult(
            read_id, f"{bc1}|{bc2}", potential_umi, s1_sc + s2_sc,
            UMI_good=umi_good, polyT=polyt_start, r1=r1_val,
        )

    def find_barcode_umi_no_polya(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        return self.find_barcode_umi(read_id, sequence)

    @staticmethod
    def result_type():
        return TenXBarcodeDetectionResult

    @classmethod
    def header(cls) -> str:
        return cls.result_type().header()


class TenXSplittingBarcodeDetector(TenXBarcodeDetector):
    MIN_REMAINING_SEQ = 50
    DEFAULT_POLYT_STEP = 100
    MIN_SPLIT_STEP = 150
    MAX_R1_POLYT_DISTANCE = 50

    def __init__(self, barcode_list: List[str]) -> None:
        super().__init__(barcode_list)
        self.tso_indexer = KmerIndexer([self.TSO], kmer_size=9)

    def _find_barcode_umi_fwd_local(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        polyt_start = find_polyt_start(sequence)
        if polyt_start == -1:
            return TenXBarcodeDetectionResult(read_id, polyT=-1)
        r1_occ = self.r1_indexer.get_occurrences(sequence[0 : polyt_start + 1])
        r1_start, r1_end = detect_exact_positions(
            sequence, 0, polyt_start + 1, self.r1_indexer.k, self.R1,
            r1_occ, min_score=11, end_delta=self.TERMINAL_MATCH_DELTA,
        )
        if r1_start is None:
            r1_start, r1_end = detect_exact_positions(
                sequence, 0, polyt_start + 1, self.r1_indexer.k, self.R1,
                r1_occ, min_score=8, end_delta=self.STRICT_TERMINAL_MATCH_DELTA,
            )
        if r1_start is None or r1_end is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)

        barcode_start = r1_end + 1
        barcode_end = r1_end + self.BARCODE_LEN_10X + 1
        potential_barcode = sequence[barcode_start : barcode_end + 1]
        matching = self.barcode_indexer.get_occurrences(
            potential_barcode, max_hits=self.max_barcodes_hits, min_kmers=self.min_matching_kmers
        )
        barcode, bc_score, _, bc_end = find_candidate_with_max_score_ssw(
            matching, potential_barcode, min_score=self.min_score, score_diff=self.score_diff
        )
        if barcode is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start, r1=r1_end)

        read_barcode_end = barcode_start + bc_end - 1
        u_start = read_barcode_end + 1
        u_end = polyt_start - 1
        if u_end - u_start <= 5:
            u_end = u_start + self.UMI_LEN - 1
        umi = sequence[u_start : u_end + 1]
        good_umi = self.UMI_LEN - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LEN + self.UMI_LEN_DELTA
        if not umi:
            return TenXBarcodeDetectionResult(read_id, barcode, BC_score=bc_score, polyT=polyt_start, r1=r1_end)
        return TenXBarcodeDetectionResult(read_id, barcode, umi, bc_score, good_umi, polyT=polyt_start, r1=r1_end)

    def _find_barcode_umi_split_fwd(self, read_id: str, sequence: str, offset: int = 0) -> TenXSplitBarcodeDetectionResult:
        base_res = self._find_barcode_umi_fwd_local(read_id, sequence)
        tso_start = -1
        if base_res.polyT != -1 and base_res.is_valid():
            ts0 = base_res.polyT
            ts1 = len(sequence) - 1
            tso_occ = self.tso_indexer.get_occurrences_substr(sequence, ts0, ts1)
            ts_s, _ = detect_first_exact_positions(
                sequence, ts0, ts1 + 1, self.tso_indexer.k, self.TSO, tso_occ, min_score=18, start_delta=3, end_delta=3
            )
            if ts_s is not None:
                tso_start = ts_s

        result = TenXSplitBarcodeDetectionResult(
            read_id, barcode=base_res.get_barcode(), UMI=base_res.get_umi(),
            BC_score=base_res.BC_score, UMI_good=base_res.UMI_good,
            polyT=base_res.polyT, r1=base_res.r1, tso=tso_start,
        )
        if offset > 0:
            result.update_coordinates(offset)
        return result

    def _is_consistent_detection(self, r: TenXSplitBarcodeDetectionResult) -> bool:
        if r.r1 == -1 or r.polyT == -1 or not r.is_valid():
            return True
        return r.r1 < r.polyT and (r.polyT - r.r1) <= self.MAX_R1_POLYT_DISTANCE

    def _scan_strand(self, read_id: str, sequence: str, strand: str, read_result: SplittingBarcodeDetectionResult) -> None:
        current_start = 0
        prev_start = -1
        while True:
            r = self._find_barcode_umi_split_fwd(read_id, sequence[current_start:], offset=current_start)
            if r.polyT == -1:
                break
            r.set_strand(strand)
            if self._is_consistent_detection(r):
                read_result.append(r)
                next_start = (r.tso + len(self.TSO)) if r.tso != -1 else (r.polyT + self.DEFAULT_POLYT_STEP)
            else:
                next_start = r.r1 - len(self.R1) - 10
            next_start = max(prev_start + self.MIN_SPLIT_STEP, next_start)
            prev_start = next_start
            current_start = next_start
            if len(sequence) - current_start < self.MIN_REMAINING_SEQ:
                break

    def find_barcode_umi(self, read_id: str, sequence: str) -> SplittingBarcodeDetectionResult:
        read_result = SplittingBarcodeDetectionResult(read_id)
        self._scan_strand(read_id, sequence, "+", read_result)
        rev_seq = reverese_complement(sequence)
        self._scan_strand(read_id, rev_seq, "-", read_result)

        if not any(r.is_valid() for r in read_result.detected_patterns):
            for try_seq, strand in ((sequence, "+"), (rev_seq, "-")):
                polyt = find_polyt_start(try_seq)
                if polyt != -1:
                    base_res = self._find_barcode_near_polyt(read_id, try_seq, polyt)
                    if base_res.is_valid():
                        sr = TenXSplitBarcodeDetectionResult(
                            read_id, base_res.get_barcode(), base_res.get_umi(),
                            base_res.BC_score, base_res.UMI_good, polyT=base_res.polyT, r1=-1, tso=-1,
                        )
                        sr.set_strand(strand)
                        read_result.detected_patterns = [sr]
                        break
        if read_result.empty():
            read_result.append(self._find_barcode_umi_split_fwd(read_id, sequence))
        read_result.filter()
        return read_result

    @staticmethod
    def result_type():
        return TenXSplitBarcodeDetectionResult

    @classmethod
    def header(cls) -> str:
        return cls.result_type().header()


class TenXv2SplittingBarcodeDetector(TenXSplittingBarcodeDetector):
    UMI_LEN = 10
    TSO = "GTACTCTGCGTTGATACCACTGCTT"
