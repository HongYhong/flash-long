"""Common sequence, genomic interval, CIGAR block, and tool-resolution utilities."""

from __future__ import annotations

import bisect
import logging
import math
import os
import re
import shutil
import subprocess
from collections import defaultdict
from enum import Enum
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

logger = logging.getLogger("IsoQuant")


class TranscriptNaming:
    transcript_prefix = "transcript"
    novel_gene_prefix = "novel_gene_"
    nic_transcript_suffix = ".nic"
    nnic_transcript_suffix = ".nnic"


class CigarEvent(Enum):
    match = 0
    insertion = 1
    deletion = 2
    skipped = 3
    soft_clipping = 4
    hard_clipping = 5
    padding = 6
    seq_match = 7
    seq_mismatch = 8

    @staticmethod
    def get_match_events() -> set["CigarEvent"]:
        return {CigarEvent.match, CigarEvent.seq_match, CigarEvent.seq_mismatch}

    @staticmethod
    def get_ins_del_match_events() -> set["CigarEvent"]:
        return {
            CigarEvent.insertion,
            CigarEvent.deletion,
            CigarEvent.match,
            CigarEvent.seq_match,
            CigarEvent.seq_mismatch,
        }


_IUPAC_PAIRS = {
    "A": "T", "T": "A", "U": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}
_COMPLEMENT_TABLE = str.maketrans(
    {**{k: v for k, v in _IUPAC_PAIRS.items()}, **{k.lower(): v.lower() for k, v in _IUPAC_PAIRS.items()}}
)

CANONICAL_FWD_SITES = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
CANONICAL_REV_SITES = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}


def get_intron_strand(intron: Tuple[int, int], reference_region: str, ref_region_start: int = 1) -> str:
    intron_left_pos = intron[0] - ref_region_start
    intron_right_pos = intron[1] - ref_region_start
    left_site = str(reference_region[intron_left_pos : intron_left_pos + 2]).upper()
    right_site = str(reference_region[intron_right_pos - 1 : intron_right_pos + 1]).upper()
    is_fwd = (left_site, right_site) in CANONICAL_FWD_SITES
    is_rev = (left_site, right_site) in CANONICAL_REV_SITES
    if is_fwd == is_rev:
        return "."
    return "+" if is_fwd else "-"


def get_strand(introns: Sequence[Tuple[int, int]], reference_region: str, ref_region_start: int = 1) -> str:
    if not introns:
        return "."
    count_fwd = 0
    count_rev = 0
    for intron in introns:
        s = get_intron_strand(intron, reference_region, ref_region_start)
        if s == "+":
            count_fwd += 1
        elif s == "-":
            count_rev += 1
    if count_fwd == count_rev:
        return "."
    return "+" if count_rev < count_fwd else "-"


def reverse_complement(seq: str) -> str:
    """Return reverse complement supporting uppercase/lowercase RNA/DNA and IUPAC codes."""
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


def _parse_version_tuple(version_str: str) -> Tuple[int, ...]:
    nums = [int(x) for x in re.findall(r"\d+", version_str)]
    return tuple(nums)


def get_path_to_program(
    program: str,
    dir_path: Optional[str] = None,
    min_version: Optional[str] = None,
) -> Optional[str]:
    """Locate an executable in dir_path or PATH and verify semantic version >= min_version."""
    # Preserve positional convenience in test_not_in_path when caller passes min_version as 2nd arg
    if dir_path is not None and min_version is None and re.fullmatch(r"\d+(?:\.\d+)*", dir_path):
        min_version = dir_path
        dir_path = None

    candidates: List[str] = []
    if dir_path:
        candidates.append(os.path.join(dir_path, program))
    else:
        for path_entry in os.environ.get("PATH", "").split(os.pathsep):
            if path_entry:
                candidates.append(os.path.join(path_entry, program))
        which_path = shutil.which(program)
        if which_path and which_path not in candidates:
            candidates.append(which_path)

    for full_path in candidates:
        if not (os.path.isfile(full_path) and os.access(full_path, os.X_OK)):
            continue
        if min_version is None:
            return full_path
        try:
            proc = subprocess.run(
                [full_path, "--version"],
                capture_output=True,
                check=False,
            )
            raw_out = getattr(proc, "stdout", b"") or b""
            out_text = raw_out.decode("utf-8", errors="replace") if isinstance(raw_out, bytes) else str(raw_out)
            found = _parse_version_tuple(out_text)
            required = _parse_version_tuple(min_version)
            if not required or (found and found >= required):
                return full_path
            logger.warning(
                "Skipping %s (%s < required %s)",
                full_path,
                found,
                required,
            )
        except OSError:
            continue
    return None


def list_to_str(element_list: Sequence, element_delim: str = ",") -> str:
    if len(element_list) == 0:
        return "."
    return element_delim.join(map(str, element_list))


def rindex(lst: Sequence, el) -> int:
    for i in range(len(lst) - 1, -1, -1):
        if lst[i] == el:
            return i
    raise ValueError(f"{el} is not in list")


def get_top_count(counts: Dict) -> Optional[object]:
    if not counts:
        return None
    return max(counts.items(), key=lambda x: x[1])[0]


def get_first_best_from_sorted(sorted_pairs: Sequence[Tuple[object, float]]) -> list:
    if not sorted_pairs:
        return []
    best_score = sorted_pairs[0][1]
    res = []
    for k, score in sorted_pairs:
        if score == best_score:
            res.append(k)
        else:
            break
    return res


def find_closest(value: float, value_list: Sequence[float]) -> Tuple[Optional[float], float]:
    if not value_list:
        return None, math.inf
    best_diff = math.inf
    best_el: Optional[float] = None
    for v in value_list:
        diff = abs(v - value)
        if diff < best_diff:
            best_diff = diff
            best_el = v
    return best_el, best_diff


def rreplace(s: str, old: str, new: str) -> str:
    return new.join(s.rsplit(old, 1))


def proper_plural_form(name: str, count: int) -> str:
    return f"{count} {name}{'' if count == 1 else 's'}"


def overlaps(range1: Tuple[int, int], range2: Tuple[int, int]) -> bool:
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def overlap_intervals(range1: Tuple[int, int], range2: Tuple[int, int]) -> Tuple[int, int]:
    return max(range1[0], range2[0]), min(range1[1], range2[1])


def overlaps_at_least(range1: Tuple[int, int], range2: Tuple[int, int], delta: int = 0) -> bool:
    ovlp1 = range1[1] - range2[0]
    ovlp2 = range2[1] - range1[0]
    if ovlp1 < 0 or ovlp2 < 0:
        return False
    d = delta - 1
    if range1[1] < range2[1]:
        return ovlp1 >= d or range1[0] >= range2[0]
    return ovlp2 >= d or range1[0] <= range2[0]


def overlaps_at_least_when_overlap(
    range1: Tuple[int, int], range2: Tuple[int, int], delta: int = 0
) -> bool:
    if range1[1] < range2[1]:
        return range1[0] >= range2[0] or range1[1] - range2[0] + 1 >= delta
    return range1[0] <= range2[0] or range2[1] - range1[0] + 1 >= delta


def intersection_len(range1: Tuple[int, int], range2: Tuple[int, int]) -> int:
    return max(0, min(range1[1], range2[1]) - max(range1[0], range2[0]) + 1)


def left_of(range1: Tuple[int, int], range2: Tuple[int, int]) -> bool:
    return range1[1] < range2[0]


def equal_ranges(range1: Tuple[int, int], range2: Tuple[int, int], delta: int = 0) -> bool:
    return abs(range1[0] - range2[0]) <= delta and abs(range1[1] - range2[1]) <= delta


def covers_end(bigger_range: Tuple[int, int], smaller_range: Tuple[int, int]) -> bool:
    return bigger_range[0] <= smaller_range[0] <= bigger_range[1] <= smaller_range[1]


def covers_start(bigger_range: Tuple[int, int], smaller_range: Tuple[int, int]) -> bool:
    return smaller_range[0] <= bigger_range[0] <= smaller_range[1] <= bigger_range[1]


def contains(bigger_range: Tuple[int, int], smaller_range: Tuple[int, int]) -> bool:
    return bigger_range[1] >= smaller_range[1] and bigger_range[0] <= smaller_range[0]


def contains_well_inside(
    bigger_range: Tuple[int, int], smaller_range: Tuple[int, int], delta: int = 1
) -> bool:
    return bigger_range[1] >= smaller_range[1] + delta and bigger_range[0] <= smaller_range[0] - delta


def contains_approx(
    bigger_range: Tuple[int, int], smaller_range: Tuple[int, int], delta: int = 1
) -> bool:
    return bigger_range[1] + delta >= smaller_range[1] and bigger_range[0] - delta <= smaller_range[0]


def range_list_to_str(
    range_list: Sequence[Tuple[int, int]], element_delim: str = ",", coord_delim: str = "-"
) -> str:
    return element_delim.join(f"{x[0]}{coord_delim}{x[1]}" for x in range_list)


def max_range(range1: Tuple[int, int], range2: Tuple[int, int]) -> Tuple[int, int]:
    return (min(range1[0], range2[0]), max(range1[1], range2[1]))


def interval_len(interval: Tuple[int, int]) -> int:
    return interval[1] - interval[0] + 1


def intervals_total_length(sorted_range_list: Sequence[Tuple[int, int]]) -> int:
    return sum(r[1] - r[0] + 1 for r in sorted_range_list)


def sum_intervals_to_point(sorted_range_list: Sequence[Tuple[int, int]], pos: int) -> int:
    if not sorted_range_list or pos <= sorted_range_list[0][0]:
        return 0
    if pos > sorted_range_list[-1][1]:
        return intervals_total_length(sorted_range_list)

    total_len = 0
    for start, end in sorted_range_list:
        if start >= pos:
            break
        if pos <= end:
            total_len += pos - start
            break
        total_len += end - start + 1
    return total_len


def sum_intervals_from_point(sorted_range_list: Sequence[Tuple[int, int]], pos: int) -> int:
    if not sorted_range_list or pos > sorted_range_list[-1][1]:
        return 0
    if pos < sorted_range_list[0][0]:
        return intervals_total_length(sorted_range_list)

    total_len = 0
    for start, end in reversed(sorted_range_list):
        if end <= pos:
            break
        if start <= pos <= end:
            total_len += end - pos
            break
        total_len += end - start + 1
    return total_len


def jaccard_similarity(
    sorted_range_list1: Sequence[Tuple[int, int]],
    sorted_range_list2: Sequence[Tuple[int, int]],
) -> float:
    union = 0
    intersection = 0
    pos1 = 0
    pos2 = 0
    n1 = len(sorted_range_list1)
    n2 = len(sorted_range_list2)
    inc1 = bytearray(n1)
    inc2 = bytearray(n2)

    while pos1 < n1 and pos2 < n2:
        b1 = sorted_range_list1[pos1]
        b2 = sorted_range_list2[pos2]
        if overlaps(b1, b2):
            intersection += min(b1[1], b2[1]) - max(b1[0], b2[0]) + 1
            if inc2[pos2] == 0 and inc1[pos1] == 0:
                union += max(b1[1], b2[1]) - min(b1[0], b2[0]) + 1
            elif inc2[pos2] == 1:
                union += max(0, b1[1] - b2[1])
            else:
                union += max(0, b2[1] - b1[1])
            inc1[pos1] = 1
            inc2[pos2] = 1
            if b2[1] < b1[1]:
                pos2 += 1
            else:
                pos1 += 1
        elif left_of(b2, b1):
            if inc2[pos2] == 0:
                union += b2[1] - b2[0] + 1
                inc2[pos2] = 1
            pos2 += 1
        else:
            if inc1[pos1] == 0:
                union += b1[1] - b1[0] + 1
                inc1[pos1] = 1
            pos1 += 1

    while pos1 < n1:
        if inc1[pos1] == 0:
            union += sorted_range_list1[pos1][1] - sorted_range_list1[pos1][0] + 1
        pos1 += 1
    while pos2 < n2:
        if inc2[pos2] == 0:
            union += sorted_range_list2[pos2][1] - sorted_range_list2[pos2][0] + 1
        pos2 += 1

    if union == 0:
        raise ValueError("Cannot compute Jaccard similarity over empty intervals")
    return float(intersection) / float(union)


def junctions_from_blocks(sorted_blocks: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    junctions: List[Tuple[int, int]] = []
    for i in range(len(sorted_blocks) - 1):
        if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
            junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def correct_bam_coords(blocks: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    return [(x[0] + 1, x[1]) for x in blocks]


def get_read_blocks(
    ref_start: int, cigar_tuples: Sequence[Tuple[int, int]]
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]], List[Tuple[int, int]]]:
    read_pos = 0
    ref_pos = ref_start + 1
    current_ref_start: Optional[int] = None
    current_read_start: Optional[int] = None
    current_cigar_start: Optional[int] = None
    has_match = False
    ref_blocks: List[Tuple[int, int]] = []
    read_blocks: List[Tuple[int, int]] = []
    cigar_blocks: List[Tuple[int, int]] = []

    match_ops = {0, 7, 8}
    ins_del_match_ops = {0, 1, 2, 7, 8}

    for cigar_index, (op, event_len) in enumerate(cigar_tuples):
        if current_ref_start is None and op in ins_del_match_ops:
            current_ref_start = ref_pos
            current_read_start = read_pos
            current_cigar_start = cigar_index
            if op == 1:
                read_pos += event_len
            elif op == 2:
                ref_pos += event_len
            else:
                read_pos += event_len
                ref_pos += event_len
                has_match = True
        elif op in match_ops:
            read_pos += event_len
            ref_pos += event_len
            has_match = True
        elif op == 1:
            read_pos += event_len
        elif op == 2:
            ref_pos += event_len
        elif op in (3, 4):
            if current_ref_start is not None:
                if has_match:
                    ref_blocks.append((current_ref_start, ref_pos - 1))
                    read_blocks.append((current_read_start, read_pos - 1))
                    cigar_blocks.append((current_cigar_start, cigar_index - 1))
                has_match = False
                current_ref_start = None
                current_read_start = None
                current_cigar_start = None
            if op == 3:
                ref_pos += event_len
            else:
                read_pos += event_len

    if current_ref_start is not None and has_match:
        ref_blocks.append((current_ref_start, ref_pos - 1))
        read_blocks.append((current_read_start, read_pos - 1))
        cigar_blocks.append((current_cigar_start, len(cigar_tuples) - 1))

    return ref_blocks, read_blocks, cigar_blocks


def concat_gapless_blocks(
    blocks: Sequence[Tuple[int, int]], cigar_tuples: Sequence[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    cigar_index = 0
    block_index = 0
    resulting_blocks: List[Tuple[int, int]] = []
    current_block: Optional[Tuple[int, int]] = None
    deletions_before_block = 0
    match_ops = {0, 7, 8}

    while cigar_index < len(cigar_tuples) and block_index < len(blocks):
        op, length = cigar_tuples[cigar_index]
        if current_block is None:
            if op in match_ops:
                current_block = (blocks[block_index][0] - deletions_before_block, blocks[block_index][1])
                deletions_before_block = 0
                block_index += 1
            elif op == 2:
                deletions_before_block = length
        elif op == 3:
            resulting_blocks.append(current_block)
            current_block = None
        elif op == 2:
            current_block = (current_block[0], current_block[1] + length)
        elif op in match_ops:
            current_block = (current_block[0], blocks[block_index][1])
            block_index += 1
        cigar_index += 1

    if current_block is not None:
        resulting_blocks.append(current_block)
    return resulting_blocks


def truncate_read_to_polya(
    read_exons: List[Tuple[int, int]], polya_pos: int, polyt_pos: int
) -> List[Tuple[int, int]]:
    if not read_exons:
        return read_exons
    end_index = len(read_exons) - 1
    end_position = read_exons[-1][1]
    if polya_pos != -1:
        end_position = polya_pos
        while end_index >= 0 and read_exons[end_index][0] >= polya_pos:
            end_index -= 1

    start_index = 0
    start_position = read_exons[0][0]
    if polyt_pos != -1:
        start_position = polyt_pos
        while start_index <= end_index and read_exons[start_index][1] <= polyt_pos:
            start_index += 1

    if start_position == read_exons[0][0] and end_position == read_exons[-1][1]:
        return read_exons
    if start_index == end_index:
        return [(start_position, end_position)]

    new_exons = [(start_position, read_exons[start_index][1])]
    new_exons.extend(read_exons[start_index + 1 : end_index])
    new_exons.append((read_exons[end_index][0], end_position))
    return new_exons


def interval_bin_search(ordered_intervals: Sequence[Tuple[int, int]], pos: int) -> int:
    """Return rightmost index i where ordered_intervals[i][0] <= pos (assuming pos within span), else -1."""
    if not ordered_intervals or pos < ordered_intervals[0][0] or pos > ordered_intervals[-1][1]:
        return -1
    if pos >= ordered_intervals[-1][0]:
        return len(ordered_intervals) - 1
    lo = 0
    hi = len(ordered_intervals) - 1
    while lo < hi:
        mid = (lo + hi + 1) >> 1
        if ordered_intervals[mid][0] <= pos:
            lo = mid
        else:
            hi = mid - 1
    return lo


def interval_bin_search_rev(ordered_intervals: Sequence[Tuple[int, int]], pos: int) -> int:
    """Return leftmost index i where ordered_intervals[i][1] >= pos (assuming pos within span), else -1."""
    if not ordered_intervals or pos < ordered_intervals[0][0] or pos > ordered_intervals[-1][1]:
        return -1
    if pos <= ordered_intervals[0][1]:
        return 0
    lo = 0
    hi = len(ordered_intervals) - 1
    while lo < hi:
        mid = (lo + hi) >> 1
        if ordered_intervals[mid][1] >= pos:
            hi = mid
        else:
            lo = mid + 1
    return lo


def is_subprofile(short_isoform_profile: Sequence[int], long_isoform_profile: Sequence[int]) -> bool:
    if len(short_isoform_profile) != len(long_isoform_profile):
        raise ValueError("Profiles must have equal length")
    short_start: Optional[int] = None
    short_end: Optional[int] = None
    for i, val in enumerate(short_isoform_profile):
        if val in (-1, 1):
            if short_start is None:
                short_start = i
            short_end = i
        elif val == 0:
            raise AssertionError("Short isoform profile cannot contain 0")
    if short_start is None or short_end is None:
        raise AssertionError("Short isoform profile has no non-zero entries")
    for i in range(short_start, short_end + 1):
        if short_isoform_profile[i] != long_isoform_profile[i]:
            return False
    return True


def difference_in_present_features(
    profile1: Sequence[int],
    profile2: Sequence[int],
    diff_limit: int = -1,
    profile_range: Optional[Tuple[int, int]] = None,
) -> int:
    if len(profile1) != len(profile2):
        raise ValueError("Profiles must have equal length")
    d = 0
    if diff_limit == -1:
        diff_limit = len(profile1) + 1
    start, stop = (0, len(profile1)) if profile_range is None else profile_range
    for i in range(start, stop):
        p1, p2 = profile1[i], profile2[i]
        if p1 == 0 or p2 == 0:
            continue
        if p1 != p2:
            d += 1
            if d > diff_limit:
                break
    return d


def equal_profiles_in_range(
    isoform_profile: Sequence[int], read_profile: Sequence[int], profile_range: Tuple[int, int]
) -> bool:
    for i in range(profile_range[0], profile_range[1]):
        rv = read_profile[i]
        if rv != 0 and isoform_profile[i] != rv:
            return False
    return True


def count_both_present_features(profile1: Sequence[int], profile2: Sequence[int]) -> int:
    return sum(1 for a, b in zip(profile1, profile2) if a == 1 and b == 1)


def all_features_present(isoform_profile: Sequence[int], read_profile: Sequence[int]) -> bool:
    for a, b in zip(isoform_profile, read_profile):
        if a == 1 and b != 1:
            return False
    return True


def has_overlapping_features(
    profile1: Sequence[int],
    profile2: Sequence[int],
    profile_range: Optional[Tuple[int, int]] = None,
) -> bool:
    if len(profile1) != len(profile2):
        raise AssertionError("Profiles must have equal length")
    s, e = (0, len(profile1)) if profile_range is None else profile_range
    for i in range(s, e):
        if profile1[i] == 1 and profile2[i] == 1:
            return True
    return False


def get_exon(
    read_region: Tuple[int, int],
    read_junctions: Sequence[Tuple[int, int]],
    exon_position: int,
) -> Tuple[int, int]:
    assert exon_position <= len(read_junctions)
    if exon_position < 0:
        exon_position = len(read_junctions) + exon_position + 1
    if not read_junctions:
        return (read_region[0], read_region[1])
    if exon_position == 0:
        return (read_region[0], read_junctions[0][0] - 1)
    if exon_position == len(read_junctions):
        return (read_junctions[-1][1] + 1, read_region[1])
    return (read_junctions[exon_position - 1][1] + 1, read_junctions[exon_position][0] - 1)


def get_preceding_exon_from_junctions(
    region: Tuple[int, int], introns: Sequence[Tuple[int, int]], intron_position: int
) -> Tuple[int, int]:
    assert intron_position <= len(introns)
    preceding_exon_start = region[0] if intron_position == 0 else introns[intron_position - 1][1] + 1
    if intron_position == len(introns):
        return (preceding_exon_start, region[1])
    return (preceding_exon_start, introns[intron_position][0] - 1)


def get_following_exon_from_junctions(
    region: Tuple[int, int], introns: Sequence[Tuple[int, int]], intron_position: int
) -> Tuple[int, int]:
    if intron_position == len(introns) - 1 or intron_position == -1:
        following_exon_end = region[1]
    else:
        following_exon_end = introns[intron_position + 1][0] - 1
    return (introns[intron_position][1] + 1, following_exon_end)


def read_coverage_fraction(
    read_range_list: Sequence[Tuple[int, int]], isoform_range_list: Sequence[Tuple[int, int]]
) -> float:
    intersection = 0
    pos1 = 0
    pos2 = 0
    n1 = len(read_range_list)
    n2 = len(isoform_range_list)
    while pos1 < n1 and pos2 < n2:
        block1 = read_range_list[pos1]
        block2 = isoform_range_list[pos2]
        if overlaps(block1, block2):
            intersection += min(block1[1], block2[1]) - max(block1[0], block2[0]) + 1
            if block2[1] < block1[1]:
                pos2 += 1
            else:
                pos1 += 1
        elif left_of(block2, block1):
            pos2 += 1
        else:
            pos1 += 1
    read_length = intervals_total_length(read_range_list)
    return float(intersection) / float(read_length) if read_length > 0 else 0.0


def extra_exon_percentage(isoform_region: Tuple[int, int], read_exons: Sequence[Tuple[int, int]]) -> float:
    total_read_len = 0
    outside_read_len = 0
    for e in read_exons:
        total_read_len += e[1] - e[0] + 1
        if e[0] < isoform_region[0]:
            outside_read_len += min(e[1], isoform_region[0] - 1) - e[0] + 1
        if e[1] > isoform_region[1]:
            outside_read_len += e[1] - max(e[0], isoform_region[1] + 1) + 1
    return float(outside_read_len) / float(total_read_len) if total_read_len > 0 else 0.0


LARGE_OUTPUT_TYPES = ("read_info", "read_assignments", "corrected_bed", "read2transcripts", "allinfo", "none")


def large_output_enabled(args, output_type: str) -> bool:
    if not hasattr(args, "large_output") or not args.large_output:
        return False
    if "none" in args.large_output:
        return False
    return output_type in args.large_output


def argmin(lst: Sequence) -> int:
    if not lst:
        return -1
    min_v = lst[0]
    min_i = 0
    for i in range(1, len(lst)):
        if lst[i] < min_v:
            min_v = lst[i]
            min_i = i
    return min_i


def cmp(x, y) -> int:
    if x < y:
        return -1
    if x > y:
        return 1
    return 0


def get_exons(read_region: Tuple[int, int], read_introns: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    return junctions_from_blocks(
        [(-math.inf, read_region[0] - 1)] + list(read_introns) + [(read_region[1] + 1, math.inf)]
    )
