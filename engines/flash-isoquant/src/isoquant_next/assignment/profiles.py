"""Compact int8 Cython profile constructors backed by contiguous arrays & nogil interval matching."""

from __future__ import annotations

from collections import defaultdict
from functools import partial
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

from isoquant_next.assignment._cy_profile_assign import (
    nonoverlapping_profile_cy,
    overlapping_profile_cy,
)
from isoquant_next.core.utils import (
    contains,
    equal_ranges,
    interval_bin_search,
    interval_bin_search_rev,
    junctions_from_blocks,
    overlaps,
    overlaps_at_least,
    overlaps_at_least_when_overlap,
)


class Int8ProfileView(Sequence[int]):
    """1-byte-per-entry int8 sequence view implementing Sequence protocol + buffer access + list equality."""

    __slots__ = ("_arr",)

    def __init__(self, arr: np.ndarray) -> None:
        self._arr = np.ascontiguousarray(arr, dtype=np.int8)

    @property
    def array(self) -> np.ndarray:
        return self._arr

    def __len__(self) -> int:
        return int(self._arr.shape[0])

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._arr[index].tolist()
        return int(self._arr[index])

    def __iter__(self):
        for v in self._arr.tolist():
            yield v

    def __contains__(self, val: object) -> bool:
        if isinstance(val, int):
            return bool(np.any(self._arr == val))
        return False

    def index(self, value: int, start: int = 0, stop: Optional[int] = None) -> int:
        n = int(self._arr.shape[0])
        s = max(0, start)
        e = n if stop is None else min(n, stop)
        sub = self._arr[s:e]
        hits = np.flatnonzero(sub == value)
        if hits.size == 0:
            raise ValueError(f"{value} is not in profile")
        return int(s + hits[0])

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Int8ProfileView):
            return bool(np.array_equal(self._arr, other._arr))
        if isinstance(other, (list, tuple)):
            return self._arr.tolist() == list(other)
        if isinstance(other, np.ndarray):
            return bool(np.array_equal(self._arr, other))
        return NotImplemented

    def __repr__(self) -> str:
        return repr(self._arr.tolist())

    def tolist(self) -> List[int]:
        return self._arr.tolist()


class MappedReadProfile:
    __slots__ = ("gene_profile", "read_profile", "read_features", "gene_profile_range")

    def __init__(
        self,
        gene_profile,
        read_profile,
        read_features,
        gene_profile_range: Optional[Tuple[int, int]] = None,
    ) -> None:
        self.gene_profile = (
            gene_profile if isinstance(gene_profile, (Int8ProfileView, list, tuple)) else Int8ProfileView(gene_profile)
        )
        self.read_profile = (
            read_profile
            if (read_profile is None or isinstance(read_profile, (Int8ProfileView, list, tuple)))
            else Int8ProfileView(read_profile)
        )
        self.read_features = read_features
        self.gene_profile_range = (
            (0, len(self.gene_profile)) if gene_profile_range is None else gene_profile_range
        )


class CombinedReadProfiles:
    __slots__ = (
        "read_intron_profile",
        "read_exon_profile",
        "read_split_exon_profile",
        "alignment",
        "polya_info",
        "cage_hits",
    )

    def __init__(
        self,
        read_intron_profile: Optional[MappedReadProfile],
        read_exon_profile: Optional[MappedReadProfile],
        read_split_exon_profile: MappedReadProfile,
        polya_info=None,
        cage_hits=-1,
        alignment=None,
    ) -> None:
        self.read_intron_profile = read_intron_profile
        self.read_exon_profile = read_exon_profile
        self.read_split_exon_profile = read_split_exon_profile
        self.alignment = alignment
        self.polya_info = polya_info
        self.cage_hits = cage_hits


def _extract_equal_delta(fn) -> Optional[int]:
    if fn is equal_ranges:
        return 0
    if isinstance(fn, partial) and fn.func is equal_ranges:
        if fn.keywords and "delta" in fn.keywords:
            return int(fn.keywords["delta"])
        return int(fn.args[0]) if fn.args else 0
    return None


def _extract_absence_spec(fn) -> Optional[Tuple[int, int]]:
    if fn is contains:
        return (0, 0)
    if fn is overlaps_at_least:
        return (1, 0)
    if isinstance(fn, partial) and fn.func is overlaps_at_least:
        d = int(fn.keywords["delta"]) if (fn.keywords and "delta" in fn.keywords) else (int(fn.args[0]) if fn.args else 0)
        return (1, d)
    return None


def _extract_nonoverlap_comp(fn) -> Optional[Tuple[int, int]]:
    if fn is overlaps:
        return (0, 0)
    if fn is overlaps_at_least_when_overlap:
        return (1, 0)
    if isinstance(fn, partial) and fn.func is overlaps_at_least_when_overlap:
        d = int(fn.keywords["delta"]) if (fn.keywords and "delta" in fn.keywords) else (int(fn.args[0]) if fn.args else 0)
        return (1, d)
    return None


def _pairs_to_int32_arrays(pairs: Sequence[Tuple[int, int]]) -> Tuple[np.ndarray, np.ndarray]:
    if not pairs:
        empty = np.empty(0, dtype=np.int32)
        return empty, empty
    arr = np.asarray(pairs, dtype=np.int32)
    return np.ascontiguousarray(arr[:, 0]), np.ascontiguousarray(arr[:, 1])


class OverlappingFeaturesProfileConstructor:
    """Fast Cython nogil profile constructor over potentially overlapping known introns/exons."""

    def __init__(
        self,
        known_features: Sequence[Tuple[int, int]],
        gene_region: Tuple[int, int],
        comparator=partial(equal_ranges, delta=0),
        absence_condition=contains,
        delta: int = 0,
    ) -> None:
        self.known_features = list(known_features)
        self.gene_region = gene_region
        self.comparator = comparator
        self.absence_condition = absence_condition
        self.delta = int(delta)
        self._ref_s, self._ref_e = _pairs_to_int32_arrays(self.known_features)
        self._comp_delta = _extract_equal_delta(comparator)
        self._abs_spec = _extract_absence_spec(absence_condition)

    @staticmethod
    def match_delta(feature1: Tuple[int, int], feature2: Tuple[int, int]) -> int:
        return abs(feature1[0] - feature2[0]) + abs(feature1[1] - feature2[1])

    def construct_intron_profile(
        self, sorted_blocks: Sequence[Tuple[int, int]], polya_position: int = -1, polyt_position: int = -1
    ) -> MappedReadProfile:
        mapped_region = (sorted_blocks[0][0], sorted_blocks[-1][1])
        read_introns = junctions_from_blocks(sorted_blocks)
        return self.construct_profile_for_features(read_introns, mapped_region, polya_position, polyt_position)

    def construct_exon_profile(
        self, sorted_blocks: Sequence[Tuple[int, int]], polya_position: int = -1, polyt_position: int = -1
    ) -> MappedReadProfile:
        mapped_region = (sorted_blocks[0][1] + self.delta, sorted_blocks[-1][0] - self.delta)
        return self.construct_profile_for_features(sorted_blocks, mapped_region, polya_position, polyt_position)

    def match_genomic_features(self, read_features: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
        gene_pos = 0
        read_pos = 0
        matched_features = defaultdict(list)
        while gene_pos < len(self.known_features) and read_pos < len(read_features):
            rf = read_features[read_pos]
            gf = self.known_features[gene_pos]
            if rf[1] < gf[0]:
                read_pos += 1
            elif gf[1] < rf[0]:
                gene_pos += 1
            elif self.comparator(rf, gf):
                matched_features[read_pos].append(gene_pos)
                gene_pos += 1
            elif overlaps(rf, gf):
                gene_pos += 1

        for rp, hits in matched_features.items():
            if len(hits) > 1:
                deltas = [self.match_delta(read_features[rp], self.known_features[gp]) for gp in hits]
                best = min(deltas)
                matched_features[rp] = [gp for gp, d in zip(hits, deltas) if d == best]

        return [
            self.known_features[matched_features[i][0]] if i in matched_features else feat
            for i, feat in enumerate(read_features)
        ]

    def construct_profile_for_features(
        self,
        read_features: Sequence[Tuple[int, int]],
        mapped_region: Tuple[int, int] = (0, 0),
        polya_position: int = -1,
        polyt_position: int = -1,
    ) -> MappedReadProfile:
        if self._comp_delta is not None and self._abs_spec is not None:
            rd_s, rd_e = _pairs_to_int32_arrays(read_features)
            gp_arr, rp_arr, _fg_idx, g_range = overlapping_profile_cy(
                self._ref_s,
                self._ref_e,
                rd_s,
                rd_e,
                self.gene_region[0],
                self.gene_region[1],
                mapped_region[0],
                mapped_region[1],
                self._comp_delta,
                self._abs_spec[0],
                self._abs_spec[1],
                polya_position,
                polyt_position,
                self.delta,
            )
            return MappedReadProfile(gp_arr, rp_arr, list(read_features), g_range)

        # Fallback Python loop for custom callable comparators/absence functions
        read_profile = [0] * len(read_features)
        intron_profile = [0] * len(self.known_features)
        matched_features = defaultdict(list)

        for i in range(len(intron_profile)):
            if self.absence_condition(mapped_region, self.known_features[i]):
                intron_profile[i] = -1
        for i in range(len(read_profile)):
            if self.absence_condition(self.gene_region, read_features[i]):
                read_profile[i] = -1

        gene_pos = 0
        read_pos = 0
        while gene_pos < len(self.known_features) and read_pos < len(read_features):
            rf = read_features[read_pos]
            gf = self.known_features[gene_pos]
            if rf[1] < gf[0]:
                if read_profile[read_pos] == 0 and gene_pos > 0:
                    read_profile[read_pos] = -1
                read_pos += 1
            elif gf[1] < rf[0]:
                if read_pos > 0:
                    intron_profile[gene_pos] = -1
                gene_pos += 1
            elif self.comparator(rf, gf):
                intron_profile[gene_pos] = 1
                read_profile[read_pos] = 1
                matched_features[read_pos].append(gene_pos)
                gene_pos += 1
            elif overlaps(rf, gf):
                if self.absence_condition(mapped_region, gf):
                    intron_profile[gene_pos] = -1
                gene_pos += 1

        for rp, hits in matched_features.items():
            if len(hits) > 1:
                deltas = [self.match_delta(read_features[rp], self.known_features[gp]) for gp in hits]
                best = min(deltas)
                for gp, d in zip(hits, deltas):
                    if d > best:
                        intron_profile[gp] = -1

        if polya_position != -1:
            for i, f in enumerate(self.known_features):
                if f[0] > polya_position + self.delta:
                    intron_profile[i] = -2
        if polyt_position != -1:
            for i, f in enumerate(self.known_features):
                if f[1] < polyt_position - self.delta:
                    intron_profile[i] = -2

        start_pos = 0
        while start_pos < len(intron_profile) and intron_profile[start_pos] == 0:
            start_pos += 1
        end_pos = len(intron_profile) - 1
        while end_pos >= 0 and intron_profile[end_pos] == 0:
            end_pos -= 1

        return MappedReadProfile(
            np.asarray(intron_profile, dtype=np.int8),
            np.asarray(read_profile, dtype=np.int8),
            list(read_features),
            (start_pos, end_pos + 1),
        )


class NonOverlappingFeaturesProfileConstructor:
    """Fast Cython nogil profile constructor over disjoint sorted split exons."""

    def __init__(
        self, known_exons: Sequence[Tuple[int, int]], comparator=overlaps, delta: int = 0
    ) -> None:
        self.known_exons = list(known_exons)
        self.comparator = comparator
        self.delta = int(delta)
        self._ref_s, self._ref_e = _pairs_to_int32_arrays(self.known_exons)
        self._comp_spec = _extract_nonoverlap_comp(comparator)

    def construct_profile(
        self, sorted_blocks: Sequence[Tuple[int, int]], polya_position: int = -1, polyt_position: int = -1
    ) -> MappedReadProfile:
        if self._comp_spec is not None:
            rd_s, rd_e = _pairs_to_int32_arrays(sorted_blocks)
            gp_arr, rp_arr, g_range = nonoverlapping_profile_cy(
                self._ref_s,
                self._ref_e,
                rd_s,
                rd_e,
                self._comp_spec[0],
                self._comp_spec[1],
                self.delta,
                polya_position,
                polyt_position,
            )
            return MappedReadProfile(gp_arr, rp_arr, list(sorted_blocks), g_range)

        exon_profile = [0] * len(self.known_exons)
        read_profile = [0] * len(sorted_blocks)
        gene_pos = 0
        read_pos = 0
        while gene_pos < len(self.known_exons) and read_pos < len(sorted_blocks):
            re = sorted_blocks[read_pos]
            ge = self.known_exons[gene_pos]
            if re[1] < ge[0]:
                if gene_pos > 0 and read_profile[read_pos] == 0:
                    read_profile[read_pos] = -1
                read_pos += 1
            elif ge[1] < re[0]:
                if read_pos > 0 and exon_profile[gene_pos] == 0:
                    exon_profile[gene_pos] = -1
                gene_pos += 1
            elif self.comparator(re, ge):
                exon_profile[gene_pos] = 1
                read_profile[read_pos] = 1
                if re[1] < ge[1]:
                    read_pos += 1
                else:
                    gene_pos += 1
            else:
                if re[1] < ge[1]:
                    read_pos += 1
                else:
                    gene_pos += 1

        if polya_position != -1:
            pidx = interval_bin_search(self.known_exons, polya_position + self.delta)
            if pidx != -1:
                for i in range(pidx + 1, len(self.known_exons)):
                    exon_profile[i] = -2
        if polyt_position != -1:
            pidx = interval_bin_search_rev(self.known_exons, polyt_position - self.delta)
            if pidx != -1:
                for i in range(0, pidx):
                    exon_profile[i] = -2

        start_pos = 0
        while start_pos < len(exon_profile) and exon_profile[start_pos] == 0:
            start_pos += 1
        end_pos = len(exon_profile) - 1
        while end_pos >= 0 and exon_profile[end_pos] == 0:
            end_pos -= 1

        return MappedReadProfile(
            np.asarray(exon_profile, dtype=np.int8),
            np.asarray(read_profile, dtype=np.int8),
            list(sorted_blocks),
            (start_pos, end_pos + 1),
        )


class CombinedProfileConstructor:
    def __init__(self, gene_info, params) -> None:
        self.gene_info = gene_info
        self.params = params
        gene_region = (gene_info.start, gene_info.end)
        self.intron_profile_constructor = OverlappingFeaturesProfileConstructor(
            self.gene_info.intron_profiles.features,
            gene_region,
            comparator=partial(equal_ranges, delta=self.params.delta),
            absence_condition=partial(overlaps_at_least, delta=self.params.minimal_intron_absence_overlap),
            delta=self.params.delta,
        )
        self.exon_profile_constructor = OverlappingFeaturesProfileConstructor(
            self.gene_info.exon_profiles.features,
            gene_region,
            comparator=partial(equal_ranges, delta=self.params.delta),
            delta=self.params.delta,
        )
        self.split_exon_profile_constructor = NonOverlappingFeaturesProfileConstructor(
            self.gene_info.split_exon_profiles.features,
            comparator=partial(overlaps_at_least_when_overlap, delta=self.params.minimal_exon_overlap),
            delta=self.params.delta,
        )

    def construct_profiles(self, sorted_blocks, polya_info, cage_hits) -> CombinedReadProfiles:
        intron_profile = self.intron_profile_constructor.construct_intron_profile(
            sorted_blocks, polya_info.external_polya_pos, polya_info.external_polyt_pos
        )
        exon_profile = (
            self.exon_profile_constructor.construct_exon_profile(
                sorted_blocks, polya_info.external_polya_pos, polya_info.external_polyt_pos
            )
            if self.params.count_exons
            else None
        )
        split_exon_profile = self.split_exon_profile_constructor.construct_profile(
            sorted_blocks, polya_info.external_polya_pos, polya_info.external_polyt_pos
        )
        return CombinedReadProfiles(
            intron_profile, exon_profile, split_exon_profile, polya_info=polya_info, cage_hits=cage_hits
        )
