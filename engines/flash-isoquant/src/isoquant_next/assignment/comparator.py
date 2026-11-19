"""Splice-junction contradiction classifier between aligned reads and candidate reference isoforms."""

from __future__ import annotations

import logging
from typing import List, Sequence, Tuple

from isoquant_next.assignment.types import (
    MatchEvent,
    MatchEventSubtype,
    SupplementaryMatchConstants,
    alternative_sites,
)
from isoquant_next.core.utils import (
    contains,
    contains_approx,
    contains_well_inside,
    equal_ranges,
    get_exon,
    get_following_exon_from_junctions,
    get_preceding_exon_from_junctions,
    interval_len,
    left_of,
    overlaps,
    overlaps_at_least,
)

logger = logging.getLogger("IsoQuant")


class JunctionComparator:
    """Classifies splicing differences (FSM/ISM/AS/IR/exon skipping/alt sites) cleanly."""

    def __init__(self, params, intron_profile_constructor) -> None:
        self.params = params
        self.intron_profile_constructor = intron_profile_constructor

    def compare_junctions(
        self,
        read_junctions: Sequence[Tuple[int, int]],
        read_region: Tuple[int, int],
        isoform_junctions: Sequence[Tuple[int, int]],
        isoform_region: Tuple[int, int],
    ) -> List[MatchEvent]:
        if not read_junctions:
            return self.get_mono_exon_subtype(read_region, isoform_junctions)

        read_pos = 0
        isoform_pos = 0
        read_present = [0] * len(read_junctions)
        iso_present = [0] * len(isoform_junctions)
        contradictory_pairs = []
        absent_region = (SupplementaryMatchConstants.absent_position, SupplementaryMatchConstants.absent_position)
        cur_region = absent_region

        while read_pos < len(read_junctions) and isoform_pos < len(isoform_junctions):
            if equal_ranges(isoform_junctions[isoform_pos], read_junctions[read_pos], self.params.delta):
                read_present[read_pos] = 1
                iso_present[isoform_pos] = 1
                if cur_region != absent_region:
                    contradictory_pairs.append(cur_region)
                    cur_region = absent_region
                read_pos += 1
                isoform_pos += 1
            elif overlaps(isoform_junctions[isoform_pos], read_junctions[read_pos]):
                read_present[read_pos] = -1
                iso_present[isoform_pos] = -1
                if cur_region == absent_region:
                    cur_region = ((read_pos, read_pos), (isoform_pos, isoform_pos))
                else:
                    cur_region = ((cur_region[0][0], read_pos), (cur_region[1][0], isoform_pos))
                if read_junctions[read_pos][1] < isoform_junctions[isoform_pos][1]:
                    read_pos += 1
                else:
                    isoform_pos += 1
            elif left_of(isoform_junctions[isoform_pos], read_junctions[read_pos]):
                if cur_region != absent_region:
                    contradictory_pairs.append(cur_region)
                    cur_region = absent_region
                if read_pos > 0 or overlaps(read_region, isoform_junctions[isoform_pos]):
                    if iso_present[isoform_pos] != -1:
                        contradictory_pairs.append(
                            ((SupplementaryMatchConstants.absent_position, read_pos), (isoform_pos, isoform_pos))
                        )
                    iso_present[isoform_pos] = -1
                isoform_pos += 1
            else:
                if cur_region != absent_region:
                    contradictory_pairs.append(cur_region)
                    cur_region = absent_region
                if isoform_pos > 0 or overlaps(isoform_region, read_junctions[read_pos]):
                    if read_present[read_pos] != -1:
                        contradictory_pairs.append(
                            ((read_pos, read_pos), (SupplementaryMatchConstants.absent_position, isoform_pos))
                        )
                    read_present[read_pos] = -1
                read_pos += 1

        if cur_region != absent_region:
            contradictory_pairs.append(cur_region)

        while read_pos < len(read_junctions):
            if overlaps(isoform_region, read_junctions[read_pos]):
                if read_present[read_pos] != -1:
                    contradictory_pairs.append(
                        ((read_pos, read_pos), (SupplementaryMatchConstants.absent_position, isoform_pos))
                    )
                    read_present[read_pos] = -1
            else:
                break
            read_pos += 1

        while isoform_pos < len(isoform_junctions):
            if overlaps(read_region, isoform_junctions[isoform_pos]):
                if iso_present[isoform_pos] != -1:
                    contradictory_pairs.append(
                        ((SupplementaryMatchConstants.absent_position, read_pos), (isoform_pos, isoform_pos))
                    )
                    iso_present[isoform_pos] = -1
            else:
                break
            isoform_pos += 1

        matching_events: List[MatchEvent] = []
        if -1 in read_present or -1 in iso_present:
            matching_events = self.detect_contradiction_type(
                read_region, read_junctions, isoform_region, isoform_junctions, contradictory_pairs
            )

        if read_present[0] == 0 or read_present[-1] == 0:
            self.add_extra_out_exon_events(matching_events, read_present, read_region, read_junctions, isoform_region[0])

        return matching_events if matching_events else [MatchEvent(MatchEventSubtype.none)]

    def detect_contradiction_type(
        self, read_region, read_junctions, isoform_region, isoform_junctions, contradictory_region_pairs
    ) -> List[MatchEvent]:
        events: List[MatchEvent] = []
        for r_c, i_c in contradictory_region_pairs:
            ev = self.compare_overlapping_contradictional_regions(
                read_region, read_junctions, isoform_region, isoform_junctions, r_c, i_c
            )
            if ev is not None:
                events.append(ev)
        return events

    def compare_overlapping_contradictional_regions(
        self, read_region, read_junctions, isoform_region, isoform_junctions, read_cregion, isoform_cregion
    ) -> Optional[MatchEvent]:
        if read_cregion[0] == SupplementaryMatchConstants.absent_position:
            overlapped = isoform_junctions[isoform_cregion[0]]
            if contains(read_region, overlapped):
                if interval_len(overlapped) <= self.params.micro_intron_length and contains_well_inside(
                    get_preceding_exon_from_junctions(read_region, read_junctions, read_cregion[1]),
                    overlapped,
                    self.params.minimal_exon_overlap,
                ):
                    return MatchEvent(MatchEventSubtype.fake_micro_intron_retention, isoform_cregion, read_cregion)
                return MatchEvent(MatchEventSubtype.intron_retention, isoform_cregion, read_cregion)
            if overlaps_at_least(read_region, overlapped, self.params.minor_exon_extension):
                sub = (
                    MatchEventSubtype.incomplete_intron_retention_left
                    if overlapped[0] <= read_region[0]
                    else MatchEventSubtype.incomplete_intron_retention_right
                )
                return MatchEvent(sub, isoform_cregion, read_cregion)
            return None

        if isoform_cregion[0] == SupplementaryMatchConstants.absent_position:
            if self.are_known_introns(read_junctions, read_cregion):
                return MatchEvent(MatchEventSubtype.extra_intron_known, isoform_cregion, read_cregion)
            if self.are_suspicious_introns(read_region, read_junctions, read_cregion):
                return MatchEvent(MatchEventSubtype.none, isoform_cregion, read_cregion)
            if (
                read_cregion[0] == 0
                and interval_len(get_exon(read_region, read_junctions, 0)) <= self.params.max_fake_terminal_exon_len
            ):
                return MatchEvent(
                    MatchEventSubtype.fake_terminal_exon_left, SupplementaryMatchConstants.extra_left_region, read_cregion
                )
            if (
                read_cregion[1] == len(read_junctions) - 1
                and interval_len(get_exon(read_region, read_junctions, -1)) <= self.params.max_fake_terminal_exon_len
            ):
                return MatchEvent(
                    MatchEventSubtype.fake_terminal_exon_right,
                    SupplementaryMatchConstants.extra_right_region,
                    read_cregion,
                )
            return MatchEvent(MatchEventSubtype.extra_intron_novel, isoform_cregion, read_cregion)

        r_tot = sum(interval_len(read_junctions[i]) for i in range(read_cregion[0], read_cregion[1] + 1))
        i_tot = sum(interval_len(isoform_junctions[i]) for i in range(isoform_cregion[0], isoform_cregion[1] + 1))
        diff = abs(r_tot - i_tot)
        sim_len = diff <= min(self.params.max_intron_abs_diff, self.params.max_intron_rel_diff * max(r_tot, i_tot))
        r_known = self.are_known_introns(read_junctions, read_cregion)

        left_ov = overlaps(
            get_exon(read_region, read_junctions, read_cregion[0]),
            get_exon(isoform_region, isoform_junctions, isoform_cregion[0]),
        )
        right_ov = overlaps(
            get_exon(read_region, read_junctions, read_cregion[1] + 1),
            get_exon(isoform_region, isoform_junctions, isoform_cregion[1] + 1),
        )
        surrounded = left_ov and right_ov

        rl, rr = read_junctions[read_cregion[0]][0], read_junctions[read_cregion[1]][1]
        il, ir = isoform_junctions[isoform_cregion[0]][0], isoform_junctions[isoform_cregion[1]][1]
        sim_l = abs(rl - il) <= 2 * self.params.delta
        sim_r = abs(rr - ir) <= 2 * self.params.delta
        sim_bounds = sim_l and sim_r
        r_in = contains_approx((il, ir), (rl, rr), self.params.delta)
        i_in = contains_approx((rl, rr), (il, ir), self.params.delta)

        event = None
        if surrounded and read_cregion[1] == read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            event = self.classify_single_intron_alternation(
                read_region, read_junctions, isoform_region, isoform_junctions, read_cregion[0], isoform_cregion[0], sim_len, r_known
            )
        elif (
            len(read_junctions) > 1
            and read_cregion[1] == read_cregion[0]
            and isoform_cregion[1] == isoform_cregion[0]
            and (
                (read_cregion[0] == 0 and isoform_cregion[0] == 0 and sim_r)
                or (
                    read_cregion[0] == len(read_junctions) - 1
                    and isoform_cregion[0] == len(isoform_junctions) - 1
                    and sim_l
                )
            )
        ):
            if read_cregion[0] == 0 and isoform_cregion[0] == 0:
                re = get_preceding_exon_from_junctions(read_region, read_junctions, 0)
                ie = get_preceding_exon_from_junctions(isoform_region, isoform_junctions, 0)
            else:
                re = get_following_exon_from_junctions(read_region, read_junctions, -1)
                ie = get_following_exon_from_junctions(isoform_region, isoform_junctions, -1)
            if abs(interval_len(re) - interval_len(ie)) < 2 * self.params.delta:
                event = (
                    MatchEventSubtype.terminal_exon_misalignment_left
                    if read_cregion[0] == 0
                    else MatchEventSubtype.terminal_exon_misalignment_right
                )
            else:
                event = (
                    MatchEventSubtype.terminal_exon_shift_known if r_known else MatchEventSubtype.terminal_exon_shift_novel
                )
        elif (
            surrounded
            and sim_bounds
            and read_cregion[1] - read_cregion[0] == isoform_cregion[1] - isoform_cregion[0] >= 1
            and diff <= 2 * self.params.delta
        ):
            event = (
                MatchEventSubtype.mutually_exclusive_exons_known
                if r_known
                else MatchEventSubtype.mutually_exclusive_exons_novel
            )
        elif surrounded and r_in and read_cregion[1] == read_cregion[0] and isoform_cregion[1] > isoform_cregion[0]:
            event = self.classify_skipped_exons(isoform_junctions, isoform_cregion, sim_len, r_known, sim_bounds)
        elif surrounded and sim_bounds and read_cregion[1] > read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            if r_known:
                event = MatchEventSubtype.exon_gain_known
            elif self.are_suspicious_introns(read_region, read_junctions, read_cregion):
                event = MatchEventSubtype.intron_retention
            else:
                event = MatchEventSubtype.exon_gain_novel
        elif surrounded and sim_len and i_in and read_cregion[1] > read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            event = MatchEventSubtype.exon_detach_known if r_known else MatchEventSubtype.exon_detach_novel

        if event is None:
            if r_known:
                event = MatchEventSubtype.alternative_structure_known
            elif surrounded and self.are_suspicious_introns(read_region, read_junctions, read_cregion):
                event = MatchEventSubtype.intron_retention
            else:
                event = MatchEventSubtype.alternative_structure_novel
        return MatchEvent(event, isoform_cregion, read_cregion)

    def classify_skipped_exons(
        self, isoform_junctions, isoform_cregion, intron_length_is_similar, read_introns_known, similar_bounds
    ):
        total_exon_len = sum(
            isoform_junctions[i + 1][0] - isoform_junctions[i][1] + 1
            for i in range(isoform_cregion[0], isoform_cregion[1])
        )
        if intron_length_is_similar:
            if total_exon_len <= self.params.max_missed_exon_len:
                return MatchEventSubtype.exon_misalignment
            return MatchEventSubtype.exon_merge_known if read_introns_known else MatchEventSubtype.exon_merge_novel
        if similar_bounds:
            return MatchEventSubtype.exon_skipping_known if read_introns_known else MatchEventSubtype.exon_skipping_novel
        return None

    def classify_single_intron_alternation(
        self, read_region, read_junctions, isoform_region, isoform_junctions, read_cpos, isoform_cpos, sim_len, r_known
    ):
        if sim_len:
            if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) <= self.params.max_intron_shift:
                return MatchEventSubtype.intron_shift
            return MatchEventSubtype.intron_migration if r_known else MatchEventSubtype.intron_alternation_novel

        event = MatchEventSubtype.intron_alternation_known if r_known else MatchEventSubtype.intron_alternation_novel
        if abs(isoform_junctions[isoform_cpos][0] - read_junctions[read_cpos][0]) <= self.params.delta:
            fe_r = get_following_exon_from_junctions(read_region, read_junctions, read_cpos)
            fe_i = get_following_exon_from_junctions(isoform_region, isoform_junctions, isoform_cpos)
            min_ov = max(1, min(self.params.min_abs_exon_overlap, round(self.params.min_rel_exon_overlap * interval_len(fe_i))))
            if overlaps_at_least(fe_r, fe_i, min_ov):
                event = alternative_sites[("right", r_known)]
        elif abs(isoform_junctions[isoform_cpos][1] - read_junctions[read_cpos][1]) <= self.params.delta:
            pe_r = get_preceding_exon_from_junctions(read_region, read_junctions, read_cpos)
            pe_i = get_preceding_exon_from_junctions(isoform_region, isoform_junctions, isoform_cpos)
            min_ov = max(1, min(self.params.min_abs_exon_overlap, round(self.params.min_rel_exon_overlap * interval_len(pe_i))))
            if overlaps_at_least(pe_r, pe_i, min_ov):
                event = alternative_sites[("left", r_known)]

        if event in {
            MatchEventSubtype.intron_alternation_novel,
            MatchEventSubtype.alt_left_site_novel,
            MatchEventSubtype.alt_right_site_novel,
        } and self.are_suspicious_introns(read_region, read_junctions, (read_cpos, read_cpos)):
            event = MatchEventSubtype.intron_retention
        return event

    def are_suspicious_introns(self, read_region, read_junctions, read_cregion) -> bool:
        tot_intron = 0
        for cpos in range(read_cregion[0], read_cregion[1] + 1):
            ilen = interval_len(read_junctions[cpos])
            if ilen > self.params.max_suspicious_intron_abs_len:
                return False
            tot_intron += ilen
        tot_exon = sum(
            interval_len(get_preceding_exon_from_junctions(read_region, read_junctions, cpos))
            for cpos in range(read_cregion[0], read_cregion[1] + 1)
        )
        tot_exon += interval_len(get_following_exon_from_junctions(read_region, read_junctions, read_cregion[-1]))
        return float(tot_intron) <= float(tot_exon) * self.params.max_suspicious_intron_rel_len

    def get_mono_exon_subtype(self, read_region, isoform_junctions) -> List[MatchEvent]:
        if not isoform_junctions:
            return [MatchEvent(MatchEventSubtype.mono_exon_match)]
        events: List[MatchEvent] = []
        for i, intron in enumerate(isoform_junctions):
            if contains(read_region, intron):
                if interval_len(intron) <= self.params.micro_intron_length and contains_well_inside(
                    read_region, intron, self.params.minimal_exon_overlap
                ):
                    events.append(
                        MatchEvent(
                            MatchEventSubtype.fake_micro_intron_retention,
                            isoform_region=(i, i),
                            read_region=(SupplementaryMatchConstants.absent_position, 0),
                        )
                    )
                else:
                    events.append(
                        MatchEvent(
                            MatchEventSubtype.unspliced_intron_retention,
                            isoform_region=(i, i),
                            read_region=(SupplementaryMatchConstants.absent_position, 0),
                        )
                    )
            elif overlaps_at_least(read_region, intron, self.params.minor_exon_extension) and not contains(intron, read_region):
                sub = (
                    MatchEventSubtype.incomplete_intron_retention_left
                    if intron[0] <= read_region[0]
                    else MatchEventSubtype.incomplete_intron_retention_right
                )
                events.append(
                    MatchEvent(sub, isoform_region=(i, i), read_region=(SupplementaryMatchConstants.absent_position, 0))
                )
        return events if events else [MatchEvent(MatchEventSubtype.mono_exonic)]

    def add_extra_out_exon_events(self, match_events, read_present, read_region, read_introns, isoform_start) -> None:
        extra_left = read_present[0] == 0
        extra_right = read_present[-1] == 0
        if all(x == 0 for x in read_present):
            if read_introns[0][0] < isoform_start:
                extra_right = False
            else:
                extra_left = False
        if extra_left:
            rp = 0
            if interval_len(get_exon(read_region, read_introns, rp)) <= self.params.max_fake_terminal_exon_len:
                match_events.append(
                    MatchEvent(MatchEventSubtype.fake_terminal_exon_left, SupplementaryMatchConstants.extra_left_region, (rp, rp))
                )
                rp += 1
            while rp < len(read_present) and read_present[rp] == 0:
                match_events.append(
                    MatchEvent(
                        MatchEventSubtype.extra_intron_flanking_left, SupplementaryMatchConstants.extra_left_region, (rp, rp)
                    )
                )
                rp += 1
        if extra_right:
            rp = len(read_present) - 1
            if interval_len(get_exon(read_region, read_introns, rp + 1)) <= self.params.max_fake_terminal_exon_len:
                match_events.append(
                    MatchEvent(
                        MatchEventSubtype.fake_terminal_exon_right, SupplementaryMatchConstants.extra_right_region, (rp, rp)
                    )
                )
                rp -= 1
            while rp >= 0 and read_present[rp] == 0:
                match_events.append(
                    MatchEvent(
                        MatchEventSubtype.extra_intron_flanking_right, SupplementaryMatchConstants.extra_right_region, (rp, rp)
                    )
                )
                rp -= 1

    def profile_for_junctions_introns(self, junctions, region):
        selected = junctions[region[0] : region[1] + 1]
        return self.intron_profile_constructor.construct_profile_for_features(selected)

    def are_known_introns(self, junctions, region) -> bool:
        prof = self.profile_for_junctions_introns(junctions, region)
        return all(el == 1 for el in prof.read_profile)
