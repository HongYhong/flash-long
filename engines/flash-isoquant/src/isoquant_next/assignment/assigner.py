"""Long-read isoform assigner using Cython int8 profile scans & exact structural classification."""

from __future__ import annotations

import logging
from collections import namedtuple
from functools import partial
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from isoquant_next.assignment._cy_profile_assign import (
    diff_present_features_cy,
    has_overlapping_features_cy,
)
from isoquant_next.assignment.comparator import JunctionComparator
from isoquant_next.assignment.profiles import Int8ProfileView, OverlappingFeaturesProfileConstructor
from isoquant_next.assignment.types import (
    IsoformMatch,
    MatchClassification,
    MatchEvent,
    MatchEventSubtype,
    PolyAVerifier,
    ReadAssignment,
    ReadAssignmentType,
    SupplementaryMatchConstants,
    elongation_cost,
    event_subtype_cost,
)
from isoquant_next.config.enums import AmbiguityResolvingMethod
from isoquant_next.core.utils import (
    contains,
    contains_approx,
    difference_in_present_features,
    equal_ranges,
    extra_exon_percentage,
    has_overlapping_features,
    jaccard_similarity,
    overlap_intervals,
    overlaps,
    read_coverage_fraction,
)

logger = logging.getLogger("IsoQuant")

IsoformDiff = namedtuple("IsoformDiff", ("isoform_id", "diff"))


def _as_int8_array(prof) -> np.ndarray:
    if isinstance(prof, Int8ProfileView):
        return prof.array
    if isinstance(prof, np.ndarray) and prof.dtype == np.int8:
        return prof
    return np.ascontiguousarray(prof, dtype=np.int8)


class LongReadAssigner:
    """Assigns spliced and unspliced long reads to annotated isoforms and categorizes events."""

    def __init__(self, gene_info, params, string_pools, quick_mode: bool = False) -> None:
        assert string_pools is not None, "string_pools is required"
        self.gene_info = gene_info
        self.params = params
        self.quick_mode = quick_mode
        self.string_pools = string_pools
        self.intron_comparator = JunctionComparator(
            params,
            OverlappingFeaturesProfileConstructor(
                self.gene_info.intron_profiles.features,
                (self.gene_info.start, self.gene_info.end),
                comparator=partial(equal_ranges, delta=self.params.delta),
            ),
        )
        self.polya_verifier = PolyAVerifier(self.gene_info, self.params)

    def get_gene_id(self, isoform_id: str) -> Optional[str]:
        if hasattr(self.gene_info, "gene_id_map"):
            return self.gene_info.gene_id_map.get(isoform_id)
        return None

    @staticmethod
    def match_profile_exact_in_range(
        read_gene_profile,
        isoform_profiles: Dict[str, Sequence[int]],
        profile_range: Optional[Tuple[int, int]],
        hint=None,
    ) -> List[IsoformDiff]:
        isoform_set = hint if hint is not None else isoform_profiles.keys()
        if len(read_gene_profile) == 0 or profile_range is None or profile_range[0] >= profile_range[1]:
            return [IsoformDiff(iid, 0) for iid in sorted(isoform_set)]
        r_arr = _as_int8_array(read_gene_profile)
        s, e = profile_range
        isoforms: List[IsoformDiff] = []
        for isoform_id in isoform_set:
            i_arr = _as_int8_array(isoform_profiles[isoform_id])
            if diff_present_features_cy(i_arr, r_arr, s, e, 0) == 0:
                isoforms.append(IsoformDiff(isoform_id, 0))
        return isoforms

    @staticmethod
    def find_matching_isoforms(
        read_profile, isoform_profiles: Dict[str, Sequence[int]], hint=None, diff_limit: int = 0
    ) -> List[str]:
        if diff_limit == 0:
            res = LongReadAssigner.match_profile_exact_in_range(
                read_profile.gene_profile,
                isoform_profiles,
                read_profile.gene_profile_range,
                hint=hint,
            )
            return sorted(x.isoform_id for x in res)
        if len(read_profile.gene_profile) == 0:
            isoform_set = hint if hint is not None else isoform_profiles.keys()
            return sorted(isoform_set)
        res = LongReadAssigner.match_profile(
            read_profile.gene_profile,
            isoform_profiles,
            hint=hint,
            diff_limit=diff_limit,
            profile_range=read_profile.gene_profile_range,
        )
        return sorted(x.isoform_id for x in res if x.diff <= diff_limit)

    @staticmethod
    def match_profile(
        read_gene_profile,
        isoform_profiles: Dict[str, Sequence[int]],
        hint=None,
        diff_limit: int = -1,
        profile_range: Optional[Tuple[int, int]] = None,
    ) -> List[IsoformDiff]:
        n_read = len(read_gene_profile)
        isoform_set = hint if hint is not None else isoform_profiles.keys()
        if n_read == 0:
            isoforms_empty: List[IsoformDiff] = []
            for isoform_id in isoform_set:
                assert len(isoform_profiles[isoform_id]) == 0, "Profile length mismatch"
                isoforms_empty.append(IsoformDiff(isoform_id, 0))
            return isoforms_empty
        r_arr = _as_int8_array(read_gene_profile)
        s, e = (0, int(r_arr.shape[0])) if profile_range is None else profile_range
        isoforms: List[IsoformDiff] = []
        for isoform_id in isoform_set:
            iso_prof = isoform_profiles[isoform_id]
            assert len(iso_prof) == n_read, "Profile length mismatch"
            i_arr = _as_int8_array(iso_prof)
            diff = int(diff_present_features_cy(i_arr, r_arr, s, e, diff_limit))
            if diff_limit == -1 or diff <= diff_limit:
                isoforms.append(IsoformDiff(isoform_id, diff))
        return sorted(isoforms, key=lambda x: x.diff)

    def find_containing_isoforms(self, read_exon_split_profile, isoform_profiles, hint=None) -> set[str]:
        isoforms: set[str] = set()
        read_exons = read_exon_split_profile.read_features
        read_region = (read_exons[0][0], read_exons[-1][1])
        isoform_set = hint if hint is not None else isoform_profiles.keys()
        for isoform_id in isoform_set:
            if contains_approx(self.gene_info.transcript_region(isoform_id), read_region, self.params.min_abs_exon_overlap):
                isoforms.add(isoform_id)
        return isoforms

    @staticmethod
    def find_overlapping_isoforms(read_exon_split_profile, isoform_profiles, hint=None) -> set[str]:
        isoforms: set[str] = set()
        r_arr = _as_int8_array(read_exon_split_profile.gene_profile)
        isoform_set = hint if hint is not None else isoform_profiles.profiles.keys()
        for isoform_id in isoform_set:
            s, e = overlap_intervals(
                read_exon_split_profile.gene_profile_range, isoform_profiles.profile_ranges[isoform_id]
            )
            if s < e:
                i_arr = _as_int8_array(isoform_profiles.profiles[isoform_id])
                if has_overlapping_features_cy(i_arr, r_arr, s, e):
                    isoforms.add(isoform_id)
        return isoforms

    def select_similar_isoforms(self, combined_read_profile):
        read_split_exon_profile = combined_read_profile.read_split_exon_profile
        overlapping = self.find_overlapping_isoforms(read_split_exon_profile, self.gene_info.split_exon_profiles)
        if not overlapping:
            return None
        sig = self.resolve_by_nucleotide_score(
            combined_read_profile,
            overlapping,
            similarity_function=self.coverage_based_nucleotide_score,
            top_scored_factor=0,
        )
        if not sig:
            return None
        intron_matching = self.match_profile(
            combined_read_profile.read_intron_profile.gene_profile,
            self.gene_info.intron_profiles.profiles,
            profile_range=combined_read_profile.read_intron_profile.gene_profile_range,
            hint=sig,
        )
        candidates = []
        read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[-1][1])
        for isoform_id, diff_introns in intron_matching:
            t_start, t_end = self.gene_info.transcript_region(isoform_id)
            extra_l = 0.5 if read_region[0] + self.params.delta < t_start else 0
            extra_r = 0.5 if read_region[1] - self.params.delta > t_end else 0
            candidates.append((isoform_id, diff_introns + extra_r + extra_l))
        best_diff = min(candidates, key=lambda x: x[1])[1]
        return [x[0] for x in candidates if x[1] <= best_diff + 3]

    def categorize_exon_elongation_subtype(self, read_split_exon_profile, isoform_id: str) -> List[MatchEvent]:
        split_exons = self.gene_info.split_exon_profiles.features
        isoform_profile = self.gene_info.split_exon_profiles.profiles[isoform_id]
        profile_range = self.gene_info.split_exon_profiles.profile_ranges[isoform_id]
        read_profile_range = read_split_exon_profile.gene_profile_range

        common_first = -1
        iso_first = profile_range[0]
        for i in range(max(iso_first, read_profile_range[0]), len(split_exons)):
            if isoform_profile[i] == 1 and read_split_exon_profile.gene_profile[i] == 1:
                common_first = i
                break

        common_last = -1
        iso_last = profile_range[1] - 1
        for i in range(min(iso_last, read_profile_range[1] - 1), -1, -1):
            if isoform_profile[i] == 1 and read_split_exon_profile.gene_profile[i] == 1:
                common_last = i
                break

        if common_first == -1 or common_last == -1:
            return []

        events: List[MatchEvent] = []
        first_exon = read_split_exon_profile.read_features[0]
        if overlaps(first_exon, split_exons[common_first]):
            extra_l = split_exons[common_first][0] - first_exon[0]
            if common_first == iso_first:
                if abs(extra_l) <= self.params.minor_exon_extension:
                    sub = (
                        MatchEventSubtype.terminal_site_match_left_precise
                        if abs(extra_l) <= self.params.delta
                        else MatchEventSubtype.terminal_site_match_left
                    )
                    events.append(MatchEvent(sub, event_info=extra_l))
                if extra_l > self.params.minor_exon_extension:
                    events.append(MatchEvent(MatchEventSubtype.major_exon_elongation_left, event_info=extra_l))
                elif extra_l > self.params.delta:
                    events.append(MatchEvent(MatchEventSubtype.exon_elongation_left, event_info=extra_l))
            elif self.params.minor_exon_extension >= extra_l > self.params.delta:
                events.append(MatchEvent(MatchEventSubtype.exon_elongation_left, event_info=extra_l))

        last_exon = read_split_exon_profile.read_features[-1]
        if overlaps(last_exon, split_exons[common_last]):
            extra_r = last_exon[1] - split_exons[common_last][1]
            if common_last == iso_last:
                if abs(extra_r) <= self.params.minor_exon_extension:
                    sub = (
                        MatchEventSubtype.terminal_site_match_right_precise
                        if abs(extra_r) <= self.params.delta
                        else MatchEventSubtype.terminal_site_match_right
                    )
                    events.append(MatchEvent(sub, event_info=extra_r))
                if extra_r > self.params.minor_exon_extension:
                    events.append(MatchEvent(MatchEventSubtype.major_exon_elongation_right, event_info=extra_r))
                elif extra_r > self.params.delta:
                    events.append(MatchEvent(MatchEventSubtype.exon_elongation_right, event_info=extra_r))
            elif self.params.minor_exon_extension >= extra_r > self.params.delta:
                events.append(MatchEvent(MatchEventSubtype.exon_elongation_right, event_info=extra_r))

        return events

    def resolve_by_nucleotide_score(
        self,
        combined_read_profile,
        matched_isoforms,
        similarity_function,
        min_similarity: float = AmbiguityResolvingMethod.minimal_score.value,
        top_scored_factor: float = AmbiguityResolvingMethod.top_scored_factor.value,
    ) -> List[str]:
        if not matched_isoforms:
            return []
        read_exons = combined_read_profile.read_split_exon_profile.read_features
        scores = [(iid, similarity_function(read_exons, self.gene_info.all_isoforms_exons[iid])) for iid in matched_isoforms]
        scores.sort(key=lambda x: x[1], reverse=True)
        best = scores[0][1]
        if top_scored_factor == 0:
            top_scored = sorted(x for x in scores if x[1] >= min_similarity)
        else:
            top_scored = sorted(x for x in scores if x[1] * top_scored_factor >= best and x[1] >= min_similarity)
        return [x[0] for x in top_scored]

    def coverage_based_nucleotide_score(self, read_exons, isoform_exons) -> float:
        cov = read_coverage_fraction(read_exons, isoform_exons)
        iso_ext = (isoform_exons[0][0] - self.params.minor_exon_extension, isoform_exons[-1][1] + self.params.minor_exon_extension)
        return cov - extra_exon_percentage(iso_ext, read_exons)

    def jaccard_based_nucleotide_score(self, read_exons, isoform_exons) -> float:
        js = jaccard_similarity(read_exons, isoform_exons)
        iso_ext = (isoform_exons[0][0] - self.params.minor_exon_extension, isoform_exons[-1][1] + self.params.minor_exon_extension)
        return js - extra_exon_percentage(iso_ext, read_exons)

    def check_read_ends(self, read_split_exon_profile, assignment: ReadAssignment) -> None:
        for match in assignment.isoform_matches:
            if match.assigned_transcript is None:
                continue
            etypes = self.categorize_exon_elongation_subtype(read_split_exon_profile, match.assigned_transcript)
            for e in etypes:
                match.add_subclassification(e)
            if any(MatchEventSubtype.is_major_elongation(e.event_type) for e in etypes):
                if not assignment.assignment_type.is_inconsistent():
                    assignment.set_assignment_type(ReadAssignmentType.inconsistent_non_intronic)
            elif any(MatchEventSubtype.is_minor_elongation(e.event_type) for e in etypes):
                if assignment.assignment_type == ReadAssignmentType.unique:
                    assignment.set_assignment_type(ReadAssignmentType.unique_minor_difference)

    def detect_ism_subtype(self, read_region: Tuple[int, int], isoform_id: str) -> MatchEvent:
        iso_introns = self.gene_info.all_isoforms_introns[isoform_id]
        iso_intron_reg = (iso_introns[0][0], iso_introns[-1][1])
        lt = iso_intron_reg[0] < read_region[0]
        rt = iso_intron_reg[1] > read_region[1]
        if lt and rt:
            et = MatchEventSubtype.ism_internal
        elif lt:
            et = MatchEventSubtype.ism_left
        elif rt:
            et = MatchEventSubtype.ism_right
        else:
            et = MatchEventSubtype.none
        return MatchEvent(et)

    def is_fsm(self, read_region: Tuple[int, int], isoform_id: str) -> bool:
        iso_introns = self.gene_info.all_isoforms_introns[isoform_id]
        return contains(read_region, (iso_introns[0][0], iso_introns[-1][1]))

    def categorize_correct_splice_match(self, combined_read_profile, isoform_id: str) -> IsoformMatch:
        read_exons = combined_read_profile.read_split_exon_profile.read_features
        read_region = (read_exons[0][0], read_exons[-1][1])
        if len(combined_read_profile.read_intron_profile.read_profile) == 0 or len(self.gene_info.all_isoforms_introns[isoform_id]) == 0:
            mc = MatchClassification.mono_exon_match
            ms = MatchEvent(MatchEventSubtype.mono_exon_match)
        elif self.is_fsm(read_region, isoform_id):
            mc = MatchClassification.full_splice_match
            ms = MatchEvent(MatchEventSubtype.fsm)
        else:
            mc = MatchClassification.incomplete_splice_match
            ms = self.detect_ism_subtype(read_region, isoform_id)
        return IsoformMatch(
            mc,
            self.string_pools,
            assigned_gene=self.get_gene_id(isoform_id),
            assigned_transcript=isoform_id,
            match_subclassification=ms,
            transcript_strand=self.gene_info.isoform_strands[isoform_id],
        )

    def categorize_multiple_splice_matches(self, combined_read_profile, isoform_ids) -> List[IsoformMatch]:
        return [self.categorize_correct_splice_match(combined_read_profile, iid) for iid in sorted(isoform_ids)]

    def categorize_correct_unspliced_match(self, combined_read_profile, isoform_id: str) -> IsoformMatch:
        iso_exons = self.gene_info.all_isoforms_exons[isoform_id]
        events = [MatchEvent(MatchEventSubtype.mono_exon_match if len(iso_exons) == 1 else MatchEventSubtype.mono_exonic)]
        mc = MatchClassification.get_mono_exon_classification(events)
        return IsoformMatch(
            mc,
            self.string_pools,
            assigned_gene=self.get_gene_id(isoform_id),
            assigned_transcript=isoform_id,
            match_subclassification=events,
            transcript_strand=self.gene_info.isoform_strands[isoform_id],
        )

    def categorize_multiple_unspliced_matches(self, combined_read_profile, isoform_ids) -> List[IsoformMatch]:
        return [self.categorize_correct_unspliced_match(combined_read_profile, iid) for iid in sorted(isoform_ids)]

    def assign_to_overlapping_genes(self, read_id: str, read_region: Tuple[int, int], match_classification) -> ReadAssignment:
        gene_matches = [
            IsoformMatch(
                match_classification,
                self.string_pools,
                assigned_gene=gid,
                transcript_strand=self.gene_info.gene_strands.get(gid, "."),
            )
            for gid, gr in self.gene_info.get_gene_regions().items()
            if overlaps(read_region, gr)
        ]
        if not gene_matches:
            gene_matches = IsoformMatch(match_classification, string_pools=self.string_pools)
        return ReadAssignment(read_id, ReadAssignmentType.noninformative, self.string_pools, match=gene_matches)

    def assign_to_isoform(self, read_id: str, combined_read_profile) -> ReadAssignment:
        read_intron_prof = combined_read_profile.read_intron_profile
        read_split_prof = combined_read_profile.read_split_exon_profile

        if self.gene_info.empty():
            return ReadAssignment(
                read_id,
                ReadAssignmentType.intergenic,
                self.string_pools,
                match=IsoformMatch(MatchClassification.intergenic, string_pools=self.string_pools),
            )

        if all(el != 1 for el in read_split_prof.read_profile) or all(
            el == 0 or el == -2 for el in read_split_prof.gene_profile
        ):
            read_region = (read_split_prof.read_features[0][0], read_split_prof.read_features[-1][1])
            gene_region = (
                self.gene_info.split_exon_profiles.features[0][0],
                self.gene_info.split_exon_profiles.features[-1][1],
            )
            if not overlaps(read_region, gene_region):
                return ReadAssignment(
                    read_id,
                    ReadAssignmentType.noninformative,
                    self.string_pools,
                    match=IsoformMatch(MatchClassification.intergenic, string_pools=self.string_pools),
                )
            if all(el != 1 for el in read_split_prof.gene_profile):
                return self.assign_to_overlapping_genes(read_id, read_region, MatchClassification.genic_intron)
            return self.assign_to_overlapping_genes(read_id, read_region, MatchClassification.genic)

        if -1 in read_intron_prof.read_profile or -1 in read_split_prof.read_profile:
            return self.match_inconsistent(read_id, combined_read_profile)
        if 0 in read_intron_prof.read_profile or 0 in read_split_prof.read_profile:
            return self.match_inconsistent(read_id, combined_read_profile)

        assignment = self.match_consistent(read_id, combined_read_profile)
        if assignment is None:
            assignment = self.match_inconsistent(read_id, combined_read_profile)
        return assignment

    def match_consistent(self, read_id: str, combined_read_profile) -> Optional[ReadAssignment]:
        read_intron_prof = combined_read_profile.read_intron_profile
        read_split_prof = combined_read_profile.read_split_exon_profile
        containing = self.find_containing_isoforms(read_split_prof, self.gene_info.split_exon_profiles.profiles)
        if not containing:
            return None
        overlapping = self.find_overlapping_isoforms(read_split_prof, self.gene_info.split_exon_profiles, hint=containing)
        if not overlapping:
            return None
        consistent = self.find_matching_isoforms(read_intron_prof, self.gene_info.intron_profiles.profiles, hint=overlapping)

        if not read_intron_prof.read_profile:
            assignment = self.match_consistent_unspliced(read_id, combined_read_profile, consistent)
        else:
            assignment = self.match_consistent_spliced(read_id, combined_read_profile, consistent)

        if assignment is not None:
            self.check_read_ends(combined_read_profile.read_split_exon_profile, assignment)
            self.verify_read_ends_for_assignment(combined_read_profile, assignment)
            if assignment.assignment_type.is_inconsistent():
                return None
        return assignment

    def match_consistent_spliced(self, read_id: str, combined_read_profile, consistent_isoforms) -> Optional[ReadAssignment]:
        iso_split = self.gene_info.split_exon_profiles.profiles
        matched = consistent_isoforms
        if len(consistent_isoforms) > 1:
            ex_matched = self.find_matching_isoforms(
                combined_read_profile.read_split_exon_profile, iso_split, hint=consistent_isoforms
            )
            if ex_matched:
                matched = ex_matched

        if len(matched) > 1:
            read_exons = combined_read_profile.read_split_exon_profile.read_features
            read_region = (read_exons[0][0], read_exons[-1][1])
            if self.params.resolve_ambiguous == AmbiguityResolvingMethod.all or (
                self.params.resolve_ambiguous == AmbiguityResolvingMethod.monoexon_and_fsm
                and any(self.is_fsm(read_region, iid) for iid in matched)
            ):
                matched = self.resolve_by_nucleotide_score(
                    combined_read_profile, matched, similarity_function=self.jaccard_based_nucleotide_score
                )

        if len(matched) == 1:
            iid = list(matched)[0]
            return ReadAssignment(
                read_id,
                ReadAssignmentType.unique,
                self.string_pools,
                match=self.categorize_correct_splice_match(combined_read_profile, iid),
            )
        if len(matched) > 1:
            return ReadAssignment(
                read_id,
                ReadAssignmentType.ambiguous,
                self.string_pools,
                match=self.categorize_multiple_splice_matches(combined_read_profile, matched),
            )
        return None

    def match_consistent_unspliced(self, read_id: str, combined_read_profile, consistent_isoforms) -> Optional[ReadAssignment]:
        matched = consistent_isoforms
        if len(consistent_isoforms) > 1 and self.params.resolve_ambiguous != AmbiguityResolvingMethod.none:
            matched = self.resolve_by_nucleotide_score(
                combined_read_profile, sorted(consistent_isoforms), similarity_function=self.jaccard_based_nucleotide_score
            )
        if len(matched) == 1:
            iid = matched[0]
            return ReadAssignment(
                read_id,
                ReadAssignmentType.unique,
                self.string_pools,
                match=self.categorize_correct_unspliced_match(combined_read_profile, iid),
            )
        if len(matched) > 1:
            return ReadAssignment(
                read_id,
                ReadAssignmentType.ambiguous,
                self.string_pools,
                match=self.categorize_multiple_unspliced_matches(combined_read_profile, matched),
            )
        return None

    def match_inconsistent(self, read_id: str, combined_read_profile) -> ReadAssignment:
        if self.quick_mode:
            return ReadAssignment(
                read_id,
                ReadAssignmentType.inconsistent,
                self.string_pools,
                match=IsoformMatch(MatchClassification.genic, string_pools=self.string_pools),
            )
        read_features = combined_read_profile.read_split_exon_profile.read_features
        read_region = (read_features[0][0], read_features[-1][1])

        best_candidates = self.select_similar_isoforms(combined_read_profile)
        if not best_candidates:
            return self.assign_to_overlapping_genes(read_id, read_region, MatchClassification.genic)
        read_matches = self.detect_inconsistensies(read_id, combined_read_profile, best_candidates)
        if not read_matches:
            return self.assign_to_overlapping_genes(read_id, read_region, MatchClassification.genic)
        best_isoforms, penalty_score = self.select_best_among_inconsistent(combined_read_profile, read_matches)
        if not best_isoforms:
            return self.assign_to_overlapping_genes(read_id, read_region, MatchClassification.genic)

        best_isoforms = sorted(best_isoforms)
        assignment_type = self.classify_assignment(best_isoforms, read_matches)
        if not combined_read_profile.read_intron_profile.read_profile:
            iso_matches = self.create_monoexon_matches(read_matches, best_isoforms)
        elif assignment_type.is_inconsistent():
            iso_matches = self.create_inconsistent_matches(read_matches, best_isoforms, penalty_score)
        else:
            iso_matches = self.create_consistent_matches(read_matches, best_isoforms, combined_read_profile)
        return ReadAssignment(read_id, assignment_type, self.string_pools, match=iso_matches)

    def classify_assignment(self, best_isoforms, read_matches) -> ReadAssignmentType:
        is_ambig = len(best_isoforms) > 1
        all_types = {e.event_type for iid in best_isoforms for e in read_matches[iid]}
        if all(MatchEventSubtype.is_consistent(e) for e in all_types):
            return ReadAssignmentType.ambiguous if is_ambig else ReadAssignmentType.unique
        if any(MatchEventSubtype.is_major_inconsistency(e) for e in all_types):
            inc_t = (
                ReadAssignmentType.inconsistent
                if any(MatchEventSubtype.is_intronic_inconsistency(e) for e in all_types)
                else ReadAssignmentType.inconsistent_non_intronic
            )
            return ReadAssignmentType.inconsistent_ambiguous if is_ambig else inc_t
        if any(MatchEventSubtype.is_minor_error(e) for e in all_types):
            return ReadAssignmentType.ambiguous if is_ambig else ReadAssignmentType.unique_minor_difference
        logger.warning("Unexpected event reported: %s", all_types)
        return ReadAssignmentType.noninformative

    def create_monoexon_matches(self, read_matches, selected_isoforms) -> List[IsoformMatch]:
        return [
            IsoformMatch(
                MatchClassification.get_mono_exon_classification(read_matches[iid]),
                self.string_pools,
                assigned_gene=self.get_gene_id(iid),
                assigned_transcript=iid,
                match_subclassification=read_matches[iid],
                transcript_strand=self.gene_info.isoform_strands[iid],
            )
            for iid in selected_isoforms
        ]

    def create_consistent_matches(self, read_matches, selected_isoforms, combined_read_profile) -> List[IsoformMatch]:
        matches = []
        for iid in selected_isoforms:
            m = self.categorize_correct_splice_match(combined_read_profile, iid)
            for e in read_matches[iid]:
                if e.event_type != MatchEventSubtype.none:
                    m.add_subclassification(e)
            matches.append(m)
        return matches

    def create_inconsistent_matches(self, read_matches, selected_isoforms, penalty_score: float = 0.0) -> List[IsoformMatch]:
        return [
            IsoformMatch(
                MatchClassification.get_inconsistency_classification(read_matches[iid]),
                self.string_pools,
                assigned_gene=self.get_gene_id(iid),
                assigned_transcript=iid,
                match_subclassification=read_matches[iid],
                transcript_strand=self.gene_info.isoform_strands[iid],
                penalty_score=penalty_score,
            )
            for iid in selected_isoforms
        ]

    def detect_inconsistensies(self, read_id: str, combined_read_profile, matched_isoforms) -> Dict[str, List[MatchEvent]]:
        read_intron_prof = combined_read_profile.read_intron_profile
        read_split_prof = combined_read_profile.read_split_exon_profile
        read_matches: Dict[str, List[MatchEvent]] = {}
        read_region = (read_split_prof.read_features[0][0], read_split_prof.read_features[-1][1])
        for iid in sorted(matched_isoforms):
            iso_introns = self.gene_info.all_isoforms_introns[iid]
            iso_region = self.gene_info.transcript_region(iid)
            events = self.intron_comparator.compare_junctions(
                read_intron_prof.read_features, read_region, iso_introns, iso_region
            )
            if len(events) == 1 and events[0].event_type == MatchEventSubtype.undefined:
                continue
            events += self.categorize_exon_elongation_subtype(read_split_prof, iid)
            events = self.polya_verifier.verify_read_ends(combined_read_profile, iid, events)
            read_matches[iid] = events
        return read_matches

    def select_best_among_inconsistent(self, combined_read_profile, read_matches) -> Tuple[List[str], float]:
        isoform_scores: List[Tuple[str, float]] = []
        for iid, match_events in read_matches.items():
            penalty = 0.0
            for e in match_events:
                cnt = 1
                if (
                    e.isoform_region != SupplementaryMatchConstants.undefined_region
                    and SupplementaryMatchConstants.absent_position not in e.isoform_region
                    and e.event_type in {MatchEventSubtype.exon_skipping_known, MatchEventSubtype.exon_skipping_novel}
                ):
                    cnt = e.isoform_region[1] - e.isoform_region[0] + 1
                elif (
                    e.read_region != SupplementaryMatchConstants.undefined_region
                    and SupplementaryMatchConstants.absent_position not in e.read_region
                ):
                    cnt = e.read_region[1] - e.read_region[0] + 1
                    if e.event_type in {
                        MatchEventSubtype.exon_gain_novel,
                        MatchEventSubtype.exon_gain_known,
                        MatchEventSubtype.mutually_exclusive_exons_novel,
                        MatchEventSubtype.mutually_exclusive_exons_known,
                        MatchEventSubtype.exon_detach_known,
                        MatchEventSubtype.exon_detach_novel,
                    }:
                        cnt -= 1
                    elif e.event_type in {
                        MatchEventSubtype.intron_retention,
                        MatchEventSubtype.unspliced_intron_retention,
                        MatchEventSubtype.fake_micro_intron_retention,
                        MatchEventSubtype.incomplete_intron_retention_left,
                        MatchEventSubtype.incomplete_intron_retention_right,
                    }:
                        cnt = 1
                cost = event_subtype_cost[e.event_type]
                if e.event_type in {
                    MatchEventSubtype.major_exon_elongation_left,
                    MatchEventSubtype.major_exon_elongation_right,
                    MatchEventSubtype.exon_elongation_right,
                    MatchEventSubtype.exon_elongation_left,
                }:
                    cost = elongation_cost(self.params, e.event_info)
                penalty += cost * cnt
            isoform_scores.append((iid, penalty))

        min_pen = min(isoform_scores, key=lambda x: x[1])[1]
        best = [iid for iid, pen in isoform_scores if pen == min_pen]
        if len(best) > 1:
            best = self.resolve_by_nucleotide_score(
                combined_read_profile, best, similarity_function=self.coverage_based_nucleotide_score
            )
        return best, min_pen

    def verify_read_ends_for_assignment(self, combined_read_profile, assignment: ReadAssignment) -> None:
        match_dict = {}
        for m in assignment.isoform_matches:
            if m.assigned_transcript is None:
                continue
            m.match_subclassifications = self.polya_verifier.verify_read_ends(
                combined_read_profile, m.assigned_transcript, m.match_subclassifications
            )
            match_dict[m.assigned_transcript] = m.match_subclassifications
        assignment.assignment_type = self.classify_assignment(match_dict.keys(), match_dict)
