"""PolyA tail and PolyT head coordinate finder sharing core/utils reverse_complement."""

from __future__ import annotations

import logging

from isoquant_next.alignment.alignment_info import PolyAInfo
from isoquant_next.core.utils import reverse_complement

logger = logging.getLogger("IsoQuant")


def move_ref_coord_alogn_alignment(alignment, shift: int) -> int:
    """Move ``shift`` query bases along CIGAR alignment from start (>0) or end (<0)."""
    if shift == 0:
        return 0

    cigar_tuples = alignment.cigartuples
    assert cigar_tuples

    direction = 1 if shift > 0 else -1
    shift = abs(shift)

    if direction == 1:
        current_pos = 0
        if len(cigar_tuples) > 1 and cigar_tuples[0][0] == 5 and cigar_tuples[1][0] == 4:
            current_pos = 2
        elif cigar_tuples[0][0] in (4, 5):
            current_pos = 1
    else:
        current_pos = -1
        if len(cigar_tuples) > 1 and cigar_tuples[-1][0] == 5 and cigar_tuples[-2][0] == 4:
            current_pos -= 2
        elif cigar_tuples[-1][0] in (4, 5):
            current_pos -= 1

    read_length_consumed = 0
    reference_length_consumed = 0
    shift += 1
    n_cigar = len(cigar_tuples)
    while -n_cigar <= current_pos < n_cigar and read_length_consumed < shift:
        cigar_event, event_len = cigar_tuples[current_pos]
        if cigar_event == 1:
            read_length_consumed += event_len
        elif cigar_event in (2, 3):
            reference_length_consumed += event_len
        elif cigar_event in (0, 7, 8):
            remaining_bases = shift - read_length_consumed
            if event_len < remaining_bases:
                reference_length_consumed += event_len
                read_length_consumed += event_len
            else:
                read_length_consumed += remaining_bases
                reference_length_consumed += remaining_bases
        elif cigar_event in (4, 5):
            break
        else:
            logger.warning("Unexpected event: %s", cigar_event)
        current_pos += direction

    return reference_length_consumed - 1


class PolyAFinder:
    """Sliding-window polyA tail / polyT head detector."""

    def __init__(self, window_size: int = 16, min_polya_fraction: float = 0.75) -> None:
        self.window_size = window_size
        self.min_polya_fraction = min_polya_fraction
        self.polyA_count = int(self.window_size * self.min_polya_fraction)

    def detect_polya(self, alignment) -> PolyAInfo:
        return PolyAInfo(
            self.find_polya_external(alignment),
            self.find_polyt_external(alignment),
            self.find_polya_internal(alignment),
            self.find_polyt_internal(alignment),
        )

    def find_polya_external(self, alignment) -> int:
        return self.find_polya_tail(alignment, 2, 2 * self.window_size)

    def find_polya_internal(self, alignment) -> int:
        return self.find_polya_tail(alignment, 4 * self.window_size, 2, check_entire_tail=True)

    def find_polyt_external(self, alignment) -> int:
        return self.find_polyt_head(alignment, 2, 2 * self.window_size)

    def find_polyt_internal(self, alignment) -> int:
        return self.find_polyt_head(alignment, 4 * self.window_size, 2, check_entire_head=True)

    def find_polya_tail(self, alignment, from_pos: int, to_pos: int, check_entire_tail: bool = False) -> int:
        cigar_tuples = alignment.cigartuples
        soft_clipped_tail_len = 0
        if len(cigar_tuples) > 1 and cigar_tuples[-1][0] == 5 and cigar_tuples[-2][0] == 4:
            soft_clipped_tail_len = cigar_tuples[-2][1]
        elif cigar_tuples[-1][0] == 4:
            soft_clipped_tail_len = cigar_tuples[-1][1]

        seq = alignment.seq
        if not seq:
            return -1
        assert soft_clipped_tail_len < len(seq)

        read_mapped_region_end = len(seq) - soft_clipped_tail_len
        to_check_start = max(0, read_mapped_region_end - from_pos)
        to_check_end = min(len(seq), read_mapped_region_end + to_pos + 1)
        sequence_to_check = seq[to_check_start:to_check_end].upper()
        pos = self.find_polya(sequence_to_check)

        if check_entire_tail and pos != -1:
            entire_tail = sequence_to_check[pos:]
            if entire_tail.count("A") < len(entire_tail) * self.min_polya_fraction:
                logger.debug("Internal polyA seems unreliable")
                pos = -1

        if pos == -1:
            logger.debug("No polyA found")
            return -1

        pos = to_check_start + pos
        if pos >= read_mapped_region_end:
            reference_polya_start = alignment.reference_end + (pos - read_mapped_region_end)
        else:
            ref_shift = move_ref_coord_alogn_alignment(alignment, pos - read_mapped_region_end)
            reference_polya_start = alignment.reference_end - ref_shift
        logger.debug("PolyA found at position %d", reference_polya_start)
        return reference_polya_start

    def find_polyt_head(self, alignment, from_pos: int, to_pos: int, check_entire_head: bool = False) -> int:
        cigar_tuples = alignment.cigartuples
        soft_clipped_head_len = 0
        if len(cigar_tuples) > 1 and cigar_tuples[0][0] == 5 and cigar_tuples[1][0] == 4:
            soft_clipped_head_len = cigar_tuples[1][1]
        elif cigar_tuples[0][0] == 4:
            soft_clipped_head_len = cigar_tuples[0][1]

        seq = alignment.seq
        if not seq:
            return -1
        assert soft_clipped_head_len < len(seq)

        read_mapped_region_start = soft_clipped_head_len
        to_check_start = max(0, read_mapped_region_start - to_pos)
        to_check_end = min(len(seq), read_mapped_region_start + from_pos + 1)
        sequence_to_check = reverse_complement(str(seq[to_check_start:to_check_end])).upper()
        pos = self.find_polya(sequence_to_check)

        if check_entire_head and pos != -1:
            entire_tail = sequence_to_check[pos:]
            if entire_tail.count("A") < len(entire_tail) * self.min_polya_fraction:
                pos = -1

        if pos == -1:
            logger.debug("No polyT found")
            return -1

        pos = to_check_end - pos - 1
        if pos <= read_mapped_region_start:
            reference_polyt_end = alignment.reference_start - (read_mapped_region_start - pos)
        else:
            ref_shift = move_ref_coord_alogn_alignment(alignment, pos - read_mapped_region_start)
            reference_polyt_end = alignment.reference_start + ref_shift
        logger.debug("PolyT found at position %d", reference_polyt_end)
        return max(1, reference_polyt_end)

    def find_polya(self, seq: str) -> int:
        n = len(seq)
        if n < self.window_size:
            return -1
        i = 0
        a_count = seq[0 : self.window_size].count("A")
        limit = n - self.window_size
        while i < limit:
            if a_count >= self.polyA_count:
                break
            first_base_a = seq[i] == "A"
            new_base_a = seq[i + self.window_size] == "A"
            if first_base_a and not new_base_a:
                a_count -= 1
            elif not first_base_a and new_base_a:
                a_count += 1
            i += 1
        if i >= limit:
            return -1
        return i + max(0, seq[i:].find("AA"))
