"""Short-read Illumina junction exon corrector using bisect lookup and context-managed BAM IO."""

from __future__ import annotations

import bisect
import logging
import math
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pysam

from isoquant_next.core.utils import junctions_from_blocks, overlaps

logger = logging.getLogger("IsoQuant")


def get_exons(
    read_region: Tuple[int, int], read_introns: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    return junctions_from_blocks(
        [(-math.inf, read_region[0] - 1)] + list(read_introns) + [(read_region[1] + 1, math.inf)]
    )


def validate_exons(novel_exons: Sequence[Tuple[int, int]]) -> bool:
    return list(novel_exons) == sorted(novel_exons) and all(0 < x[0] <= x[1] for x in novel_exons)


class IlluminaExonCorrector:
    """Corrects long-read junction boundary shifts and short skipped exons using Illumina splice junctions."""

    MAX_SCORE = 10**12
    ABSENT_INTRON = (0, 0)
    EXON_LENGTH = 50
    SIDE_DIFF = 25

    def __init__(
        self,
        chromosome: str | int,
        start: int,
        end: int,
        short_read_files: Sequence[str],
    ) -> None:
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.short_introns, self.counts = self.get_introns(short_read_files, str(chromosome), start, end)
        self._sorted_introns: List[Tuple[int, int]] = sorted(self.short_introns)
        self._sorted_starts: List[int] = [x[0] for x in self._sorted_introns]
        self._max_short_span: int = max((b - a for a, b in self._sorted_introns), default=0)

    @classmethod
    def from_data(cls, short_introns: Set[Tuple[int, int]] | Iterable[Tuple[int, int]]) -> "IlluminaExonCorrector":
        corrector = cls.__new__(cls)
        corrector.chromosome = 0
        corrector.start = 0
        corrector.end = 0
        corrector.short_introns = set(short_introns)
        corrector.counts = {}
        corrector._sorted_introns = sorted(corrector.short_introns)
        corrector._sorted_starts = [x[0] for x in corrector._sorted_introns]
        corrector._max_short_span = max((b - a for a, b in corrector._sorted_introns), default=0)
        return corrector

    def get_introns(
        self, files: Sequence[str], chromosome: str, start: int, end: int
    ) -> Tuple[Set[Tuple[int, int]], Dict[Tuple[int, int], int]]:
        introns: Dict[Tuple[int, int], int] = {}
        for f in files:
            with pysam.AlignmentFile(f, "rb") as samfile:
                intr = samfile.find_introns(samfile.fetch(chromosome, start=start, stop=end))
                self.merge_dictionaries(introns, intr)
        # Convert BAM 0-based start -> 1-based GTF intron coordinates; ignore integer sentinel keys cleanly
        i_list: Set[Tuple[int, int]] = {
            (k[0] + 1, k[1]) for k in introns.keys() if not isinstance(k, int)
        }
        return i_list, introns

    @staticmethod
    def merge_dictionaries(old: Dict, new: Dict) -> Dict:
        for entry, val in new.items():
            old[entry] = old.get(entry, 0) + val
        return old

    def correct_read(self, alignment_info) -> List[Tuple[int, int]]:
        return self.correct_exons(alignment_info.read_exons)

    @staticmethod
    def skipped_score(left: Tuple[int, int], right: Tuple[int, int], old: Tuple[int, int]) -> int:
        return (old[0] - left[0]) + (right[1] - old[1]) - (right[0] - left[1])

    @staticmethod
    def better_skipped(left: Tuple[int, int], right: Tuple[int, int], old: Tuple[int, int], score: int) -> bool:
        return IlluminaExonCorrector.skipped_score(left, right, old) < score

    @staticmethod
    def right_length(left: Tuple[int, int], right: Tuple[int, int], old: Tuple[int, int]) -> bool:
        middle = abs(right[0] - left[1]) <= IlluminaExonCorrector.EXON_LENGTH
        left_side = abs(old[0] - left[0]) <= IlluminaExonCorrector.SIDE_DIFF
        right_side = abs(right[1] - old[1]) <= IlluminaExonCorrector.SIDE_DIFF
        return middle and left_side and right_side

    @staticmethod
    def one_differs(left: Tuple[int, int], right: Tuple[int, int], old: Tuple[int, int]) -> bool:
        return left[0] != old[0] or right[1] != old[1]

    def _overlapping_short_introns(self, interval: Tuple[int, int]) -> List[Tuple[int, int]]:
        if not self._sorted_introns:
            return []
        lo = bisect.bisect_left(self._sorted_starts, interval[0] - self._max_short_span)
        hi = bisect.bisect_right(self._sorted_starts, interval[1])
        return [s for s in self._sorted_introns[lo:hi] if overlaps(interval, s)]

    def correct_exons(self, exons: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
        if not exons:
            return []
        if not self._sorted_introns:
            return list(exons)
        introns = junctions_from_blocks(exons)
        corrected_introns: List[Tuple[int, int]] = []

        for i in introns:
            score = self.MAX_SCORE
            sh = self.ABSENT_INTRON
            overlapping = self._overlapping_short_introns(i)
            for s in overlapping:
                x = abs(i[0] - s[0]) + abs(i[1] - s[1])
                if x < score:
                    score = x
                    sh = s

            appended = False
            if (i[0] == sh[0] and i[1] == sh[1] - 4) or (i[1] == sh[1] and sh[0] == i[0] - 4):
                corrected_introns.append(sh)
                appended = True

            if len(overlapping) > 1 and not appended:
                best_skip = self.MAX_SCORE
                left = self.ABSENT_INTRON
                right = self.ABSENT_INTRON
                for k in range(len(overlapping) - 1):
                    x = overlapping[k]
                    for j in range(k, len(overlapping)):
                        y = overlapping[j]
                        if x[1] < y[0]:
                            l_c, r_c = x, y
                        elif x[0] > y[1]:
                            l_c, r_c = y, x
                        else:
                            continue
                        if self.right_length(l_c, r_c, i) and self.one_differs(l_c, r_c, i):
                            sc = self.skipped_score(l_c, r_c, i)
                            if sc < best_skip:
                                left, right, best_skip = l_c, r_c, sc
                if left != self.ABSENT_INTRON:
                    corrected_introns.append(left)
                    corrected_introns.append(right)
                    appended = True

            if not appended:
                corrected_introns.append(i)

        res = get_exons((exons[0][0], exons[-1][1]), corrected_introns)
        if not validate_exons(res):
            logger.debug("Invalid exon order after Illumina correction: old=%s new=%s", introns, corrected_introns)
        return res


class VoidExonCorrector:
    """No-op Illumina exon corrector used when no short-read BAM is supplied."""

    def correct_read(self, alignment_info) -> List[Tuple[int, int]]:
        return alignment_info.read_exons

    def correct_exons(self, exons: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
        return list(exons)
