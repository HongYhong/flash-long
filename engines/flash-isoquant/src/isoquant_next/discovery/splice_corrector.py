"""Splice-site error correction for constructed novel transcript models."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from isoquant_next.core.utils import CANONICAL_FWD_SITES, CANONICAL_REV_SITES

logger = logging.getLogger("IsoQuant")


@dataclass
class SpliceSiteCase:
    """Per-splice-site deletion histogram & resolved canonical shift."""

    location_is_end: bool
    deletions: Dict[int, int] = field(default_factory=dict)
    most_common_del: int = -1


ACCEPTED_DEL_CASES: Tuple[int, ...] = (3, 4, 5, 6)
WINDOW_SIZE: int = 8
MIN_N_OF_ALIGNED_READS: int = 5
SUPPORTED_STRANDS: Tuple[str, ...] = ("+", "-")
BOUNDARY_MARGIN: int = WINDOW_SIZE + max(ACCEPTED_DEL_CASES)

_BAM_CDEL: int = 2
_REF_CONSUMING_OPS: frozenset = frozenset((0, 2, 3, 7, 8))

_FWD_LEFT_SITES: frozenset = frozenset(donor for donor, _ in CANONICAL_FWD_SITES)
_FWD_RIGHT_SITES: frozenset = frozenset(acceptor for _, acceptor in CANONICAL_FWD_SITES)
_REV_LEFT_SITES: frozenset = frozenset(donor for donor, _ in CANONICAL_REV_SITES)
_REV_RIGHT_SITES: frozenset = frozenset(acceptor for _, acceptor in CANONICAL_REV_SITES)


def extract_relevant_deletions(
    cigartuples: List[Tuple[int, int]],
    ref_start: int,
    exons: List[Tuple[int, int]],
    margin: int = BOUNDARY_MARGIN,
) -> List[Tuple[int, int]]:
    if not cigartuples or not exons:
        return []
    boundaries = [b for exon in exons for b in exon]
    deletions: List[Tuple[int, int]] = []
    ref_pos = ref_start
    for op, length in cigartuples:
        if op == _BAM_CDEL:
            del_end = ref_pos + length - 1
            if any(ref_pos - margin <= b <= del_end + margin for b in boundaries):
                deletions.append((ref_pos, length))
            ref_pos += length
        elif op in _REF_CONSUMING_OPS:
            ref_pos += length
    return deletions


def extract_splice_site_locations_within_aligned_read(
    read_start: int, read_end: int, exons: List[Tuple[int, int]]
) -> List[Tuple[int, bool]]:
    matching_locations: List[Tuple[int, bool]] = []
    for exon_start, exon_end in exons:
        if read_start <= exon_start <= read_end:
            matching_locations.append((exon_start, False))
        if read_start <= exon_end <= read_end:
            matching_locations.append((exon_end, True))
        if read_end <= exon_end:
            break
    return matching_locations


def _deleted_bases_in_window(deletions: List[Tuple[int, int]], low: int, high: int) -> int:
    total = 0
    for del_start, del_len in deletions:
        del_end = del_start + del_len - 1
        overlap = min(high, del_end) - max(low, del_start) + 1
        if overlap > 0:
            total += overlap
    return total


def count_deletions_for_splice_site_locations(
    read_start: int,
    read_end: int,
    deletions: List[Tuple[int, int]],
    exons: List[Tuple[int, int]],
    splice_site_cases: Dict[int, SpliceSiteCase],
    window_size: int = WINDOW_SIZE,
) -> None:
    for splice_site_location, location_is_end in extract_splice_site_locations_within_aligned_read(
        read_start, read_end, exons
    ):
        if location_is_end:
            low, high = splice_site_location - window_size + 1, splice_site_location
        else:
            low, high = splice_site_location, splice_site_location + window_size - 1
        count_of_deletions = _deleted_bases_in_window(deletions, low, high)
        case = splice_site_cases.get(splice_site_location)
        if case is None:
            case = SpliceSiteCase(location_is_end=location_is_end)
            splice_site_cases[splice_site_location] = case
        case.deletions[count_of_deletions] = case.deletions.get(count_of_deletions, 0) + 1


def compute_most_common_case_of_deletions(deletions: Dict[int, int], location_is_end: bool) -> int:
    max_count = max(deletions.values())
    modes = [k for k, v in deletions.items() if v == max_count]
    if len(modes) != 1:
        return -1
    return -modes[0] if location_is_end else modes[0]


def _corrected_site_is_canonical(
    chr_record, boundary: int, most_common_del: int, location_is_end: bool, strand: str
) -> bool:
    corrected = boundary + most_common_del
    if location_is_end:
        start = corrected
        valid = _FWD_LEFT_SITES if strand == "+" else _REV_LEFT_SITES
    else:
        start = corrected - 3
        valid = _FWD_RIGHT_SITES if strand == "+" else _REV_RIGHT_SITES
    if start < 0:
        return False
    try:
        site = str(chr_record[start : start + 2]).upper()
    except (KeyError, IndexError):
        return False
    return site in valid


def correct_splice_site_errors(
    splice_site_cases: Dict[int, SpliceSiteCase],
    chr_record,
    strand: str,
    min_n_of_aligned_reads: int = MIN_N_OF_ALIGNED_READS,
    accepted_del_cases: Tuple[int, ...] = ACCEPTED_DEL_CASES,
) -> List[int]:
    locations_with_errors: List[int] = []
    for splice_site_location, splice_site_data in splice_site_cases.items():
        if sum(splice_site_data.deletions.values()) < min_n_of_aligned_reads:
            continue
        most_common_del = compute_most_common_case_of_deletions(
            splice_site_data.deletions, splice_site_data.location_is_end
        )
        splice_site_data.most_common_del = most_common_del
        if abs(most_common_del) not in accepted_del_cases:
            continue
        if _corrected_site_is_canonical(
            chr_record, splice_site_location, most_common_del, splice_site_data.location_is_end, strand
        ):
            locations_with_errors.append(splice_site_location)
    return locations_with_errors


def generate_updated_exon_list(
    splice_site_cases: Dict[int, SpliceSiteCase],
    locations_with_errors: List[int],
    exons: List[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    corrections = set(locations_with_errors)
    updated_exons: List[Tuple[int, int]] = []
    for exon_start, exon_end in exons:
        if exon_start in corrections:
            exon_start += splice_site_cases[exon_start].most_common_del
        if exon_end in corrections:
            exon_end += splice_site_cases[exon_end].most_common_del
        updated_exons.append((exon_start, exon_end))
    return updated_exons
