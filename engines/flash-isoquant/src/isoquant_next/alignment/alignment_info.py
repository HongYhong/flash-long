"""Memory-compact AlignmentInfo & SoA locus container backed by Cython CIGAR block kernels."""

from __future__ import annotations

import logging
from typing import List, Optional, Sequence, Tuple

import numpy as np

from isoquant_next.alignment._cy_cigar_blocks import count_errors_in_arrays_cy, extract_cigar_blocks_cy
from isoquant_next.core.utils import interval_len

logger = logging.getLogger("IsoQuant")


def shift_polya(read_exons: Sequence[Tuple[int, int]], exon_count: int, polya_pos: int) -> int:
    if exon_count == 0 or exon_count == len(read_exons) or polya_pos == -1:
        return polya_pos
    dist_to_polya = 0
    for i in range(exon_count):
        exon = read_exons[-i - 1]
        if exon[0] > polya_pos:
            continue
        if dist_to_polya == 0:
            dist_to_polya += polya_pos - exon[0]
        else:
            dist_to_polya += interval_len(exon)
    return read_exons[-exon_count - 1][1] + dist_to_polya


def shift_polyt(read_exons: Sequence[Tuple[int, int]], exon_count: int, polyt_pos: int) -> int:
    if exon_count == 0 or exon_count == len(read_exons) or polyt_pos == -1:
        return polyt_pos
    dist_to_polya = 0
    for i in range(exon_count):
        exon = read_exons[i]
        if exon[1] < polyt_pos:
            continue
        if dist_to_polya == 0:
            dist_to_polya += exon[1] - polyt_pos
        else:
            dist_to_polya += interval_len(exon)
    return read_exons[exon_count][0] - dist_to_polya


class PolyAInfo:
    """Compact slots representation of detected internal/external polyA & polyT positions."""

    __slots__ = ("external_polya_pos", "external_polyt_pos", "internal_polya_pos", "internal_polyt_pos")

    def __init__(
        self,
        external_polya_pos: int = -1,
        external_polyt_pos: int = -1,
        internal_polya_pos: int = -1,
        internal_polyt_pos: int = -1,
    ) -> None:
        self.external_polya_pos = external_polya_pos
        self.external_polyt_pos = external_polyt_pos
        self.internal_polya_pos = internal_polya_pos
        self.internal_polyt_pos = internal_polyt_pos


class AlignmentInfo:
    """Per-alignment block, splice, and polyA holder using __slots__ and Cython CIGAR scanning."""

    REGION_TO_CHECK_LEN = 12

    __slots__ = (
        "alignment",
        "read_exons",
        "read_blocks",
        "cigar_blocks",
        "aligned_pairs",
        "aligned_pairs_start_index",
        "aligned_pairs_end_index",
        "read_start",
        "read_end",
        "polya_info",
        "exons_changed",
        "cage_hits",
        "combined_profile",
    )

    def __init__(self, alignment) -> None:
        self.alignment = alignment
        ct = alignment.cigartuples
        if ct:
            arr = np.asarray(ct, dtype=np.int32)
            ops = np.ascontiguousarray(arr[:, 0])
            lens = np.ascontiguousarray(arr[:, 1])
            self.read_exons, self.read_blocks, self.cigar_blocks = extract_cigar_blocks_cy(
                ops, lens, int(alignment.reference_start)
            )
        else:
            self.read_exons, self.read_blocks, self.cigar_blocks = [], [], []

        self.aligned_pairs = None
        self.aligned_pairs_start_index = None
        self.aligned_pairs_end_index = None
        self.read_start: int = self.read_exons[0][0] if self.read_exons else 0
        self.read_end: int = self.read_exons[-1][1] if self.read_exons else 0
        self.polya_info = PolyAInfo(-1, -1, -1, -1)
        self.exons_changed = False
        self.cage_hits: List = []
        self.combined_profile = None

    def set_aligned_pairs(self) -> None:
        self.aligned_pairs = self.alignment.get_aligned_pairs()
        self.aligned_pairs_start_index = []
        self.aligned_pairs_end_index = []
        if not self.read_exons:
            self.aligned_pairs_start_index.append(len(self.aligned_pairs))
            self.aligned_pairs_start_index.append(len(self.aligned_pairs))
            return
        exon_index = 0
        num_exons = len(self.read_exons)
        for i, (_read_pos, ref_pos) in enumerate(self.aligned_pairs):
            if ref_pos == self.read_exons[exon_index][0]:
                self.aligned_pairs_start_index.append(i)
            if ref_pos == self.read_exons[exon_index][1]:
                self.aligned_pairs_end_index.append(i)
                exon_index += 1
            if exon_index == num_exons:
                break
        self.aligned_pairs_start_index.append(len(self.aligned_pairs))
        self.aligned_pairs_start_index.append(len(self.aligned_pairs))

    def get_error_count(
        self,
        ref_start: int,
        ref_end: int,
        intron_index: Optional[int] = None,
        left_site: Optional[bool] = True,
        chr_record=None,
    ) -> Tuple[int, int]:
        if not self.aligned_pairs:
            self.set_aligned_pairs()

        if intron_index is None:
            selected_pairs = self.aligned_pairs
        elif left_site:
            idx = self.aligned_pairs_end_index[intron_index]
            selected_pairs = self.aligned_pairs[
                max(0, idx - self.REGION_TO_CHECK_LEN) : idx + self.REGION_TO_CHECK_LEN
            ]
        else:
            idx = self.aligned_pairs_start_index[intron_index + 1]
            selected_pairs = self.aligned_pairs[
                max(0, idx - self.REGION_TO_CHECK_LEN) : idx + self.REGION_TO_CHECK_LEN
            ]

        n = len(selected_pairs)
        if n == 0:
            return 0, 0
        q_arr = np.empty(n, dtype=np.int32)
        r_arr = np.empty(n, dtype=np.int64)
        for k, (qp, rp) in enumerate(selected_pairs):
            q_arr[k] = -1 if qp is None else qp
            r_arr[k] = -1 if rp is None else rp

        q_seq_str = getattr(self.alignment, "query_sequence", None)
        q_bytes = (
            np.frombuffer(q_seq_str.encode("ascii"), dtype=np.uint8)
            if (chr_record is not None and q_seq_str)
            else None
        )
        r_bytes = None
        ref_offset = 0
        if chr_record is not None and q_bytes is not None:
            # Slice only the bounded genomic span [ref_start-1, ref_end] so pyfaidx/FastaFile avoids whole-chrom copies
            valid_rp = [int(rp) for (_qp, rp) in selected_pairs if rp is not None]
            if valid_rp:
                r0 = max(0, min(valid_rp))
                r1 = max(valid_rp) + 1
                sub_ref = str(chr_record[r0:r1])
                ref_offset = r0
                r_bytes = np.frombuffer(sub_ref.encode("ascii"), dtype=np.uint8)

        return count_errors_in_arrays_cy(q_arr, r_arr, ref_start, ref_end, q_bytes, r_bytes, ref_offset)

    def add_polya_info(self, polya_finder, polya_fixer) -> None:
        self.polya_info = polya_finder.detect_polya(self.alignment)
        polya_exon_count, polyt_exon_count = polya_fixer.correct_read_info(self.read_exons, self.polya_info)

        if polya_exon_count > 0:
            self.polya_info.internal_polya_pos = shift_polya(
                self.read_exons, polya_exon_count, self.polya_info.internal_polya_pos
            )
            self.polya_info.external_polya_pos = shift_polya(
                self.read_exons, polya_exon_count, self.polya_info.external_polya_pos
            )
            self.read_exons = self.read_exons[:-polya_exon_count]
            self.read_blocks = self.read_blocks[:-polya_exon_count]
            self.cigar_blocks = self.cigar_blocks[:-polya_exon_count]
            self.exons_changed = True

        if polyt_exon_count > 0:
            self.polya_info.internal_polyt_pos = shift_polyt(
                self.read_exons, polyt_exon_count, self.polya_info.internal_polyt_pos
            )
            self.polya_info.external_polyt_pos = shift_polyt(
                self.read_exons, polyt_exon_count, self.polya_info.external_polyt_pos
            )
            self.read_exons = self.read_exons[polyt_exon_count:]
            self.read_blocks = self.read_blocks[polyt_exon_count:]
            self.cigar_blocks = self.cigar_blocks[polyt_exon_count:]
            self.exons_changed = True

        if self.exons_changed and self.read_exons:
            self.read_start = self.read_exons[0][0]
            self.read_end = self.read_exons[-1][1]

    def add_cage_info(self, cage_finder) -> None:
        self.cage_hits = cage_finder.find_cage_peak(self.alignment)

    def construct_profiles(self, profile_constructor) -> None:
        self.combined_profile = profile_constructor.construct_profiles(
            self.read_exons, self.polya_info, self.cage_hits
        )

    def get_seq(self, ref_start: int, ref_end: int, intron_region=None) -> str:
        raise AssertionError("get_seq is deprecated in IsoQuant AlignmentInfo")


class LocusAlignmentSoABuffer:
    """Compact Structure-of-Arrays read exon & polyA buffer for O(Locus) streaming processing."""

    __slots__ = (
        "read_ids",
        "read_starts",
        "read_ends",
        "strands",
        "polya_pos",
        "polyt_pos",
        "exon_ptr",
        "exon_starts",
        "exon_ends",
    )

    def __init__(self, alignments: Sequence[AlignmentInfo]) -> None:
        n = len(alignments)
        self.read_ids: Tuple[str, ...] = tuple(
            getattr(a.alignment, "query_name", "") or "" for a in alignments
        )
        self.read_starts = np.empty(n, dtype=np.int32)
        self.read_ends = np.empty(n, dtype=np.int32)
        self.strands = np.empty(n, dtype=np.int8)
        self.polya_pos = np.empty(n, dtype=np.int32)
        self.polyt_pos = np.empty(n, dtype=np.int32)
        self.exon_ptr = np.empty(n + 1, dtype=np.int64)
        self.exon_ptr[0] = 0

        ex_s: List[int] = []
        ex_e: List[int] = []
        for idx, ai in enumerate(alignments):
            self.read_starts[idx] = ai.read_start
            self.read_ends[idx] = ai.read_end
            is_rev = bool(getattr(ai.alignment, "is_reverse", False))
            self.strands[idx] = -1 if is_rev else 1
            pi = ai.polya_info
            pa = pi.external_polya_pos if pi.external_polya_pos != -1 else pi.internal_polya_pos
            pt = pi.external_polyt_pos if pi.external_polyt_pos != -1 else pi.internal_polyt_pos
            self.polya_pos[idx] = pa
            self.polyt_pos[idx] = pt
            for s, e in ai.read_exons:
                ex_s.append(s)
                ex_e.append(e)
            self.exon_ptr[idx + 1] = len(ex_s)

        self.exon_starts = np.asarray(ex_s, dtype=np.int32)
        self.exon_ends = np.asarray(ex_e, dtype=np.int32)

    def __len__(self) -> int:
        return int(self.read_starts.shape[0])
