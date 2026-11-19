"""Zero-pybedtools CAGE peak lookup using sorted per-chromosome bisect + prefix max_end."""

from __future__ import annotations

import logging
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, NamedTuple, Tuple, Union

logger = logging.getLogger("IsoQuant")


class CageHit(NamedTuple):
    """Lightweight BED interval record matching pybedtools.Interval overlap attributes."""

    chrom: str
    start: int
    end: int
    name: str = "."
    score: str = "."
    strand: str = "."


class CagePeakFinder:
    """In-memory BED interval index replacing pybedtools.BedTool with zero subprocess/C++ overhead."""

    def __init__(self, cage_file: Union[str, Path], shift_size: int = 50, window_size: int = 5) -> None:
        self.shift_size = shift_size
        self.window_size = window_size
        self._chrom_intervals: Dict[str, List[CageHit]] = {}
        self._chrom_starts: Dict[str, List[int]] = {}
        self._chrom_prefix_max_end: Dict[str, List[int]] = {}
        self.cage_peaks = self._load_cage_peaks(cage_file)

    def _load_cage_peaks(self, cage_file: Union[str, Path]) -> Dict[str, List[CageHit]]:
        by_chrom: Dict[str, List[CageHit]] = defaultdict(list)
        with open(cage_file, "r", encoding="utf-8") as fh:
            for line in fh:
                if not line or line[0] == "#" or line.startswith(("track", "browser")):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 3:
                    continue
                hit = CageHit(
                    chrom=parts[0],
                    start=int(parts[1]),
                    end=int(parts[2]),
                    name=parts[3] if len(parts) > 3 else ".",
                    score=parts[4] if len(parts) > 4 else ".",
                    strand=parts[5] if len(parts) > 5 else ".",
                )
                by_chrom[hit.chrom].append(hit)

        for chrom, hits in by_chrom.items():
            hits.sort(key=lambda h: (h.start, h.end))
            self._chrom_intervals[chrom] = hits
            self._chrom_starts[chrom] = [h.start for h in hits]
            pref: List[int] = []
            cur_max = -1
            for h in hits:
                if h.end > cur_max:
                    cur_max = h.end
                pref.append(cur_max)
            self._chrom_prefix_max_end[chrom] = pref

        return self._chrom_intervals

    def _get_search_region(self, alignment, extended: bool = False) -> Tuple[str, int, int, str]:
        contig = alignment.reference_name
        search_size = self.shift_size if extended else self.window_size
        if alignment.is_reverse:
            strand = "-"
            start = max(alignment.reference_end - self.window_size, 0)
            end = alignment.reference_end + search_size
        else:
            strand = "."
            start = max(alignment.reference_start - search_size, 0)
            end = alignment.reference_start + self.window_size
        return contig, start, end, strand

    def find_cage_peak(self, alignment) -> List[CageHit]:
        logger.debug("Searching for cage peak for %s ", alignment.query_name)
        contig, q_start, q_end, _ = self._get_search_region(alignment, extended=True)
        hits = self._chrom_intervals.get(contig)
        if not hits:
            logger.debug("No CAGE peaks found")
            return []

        starts = self._chrom_starts[contig]
        pref_max = self._chrom_prefix_max_end[contig]
        # pybedtools BedTool.all_hits treats book-ended BED intervals (start <= q_end and end >= q_start) as overlapping.
        idx = bisect_right(starts, q_end) - 1
        found: List[CageHit] = []
        while idx >= 0 and pref_max[idx] >= q_start:
            iv = hits[idx]
            if iv.end >= q_start and iv.start <= q_end:
                found.append(iv)
            idx -= 1

        if found:
            found.reverse()
            logger.debug("CAGE peaks found: %s", found)
        else:
            logger.debug("No CAGE peaks found")
        return found
