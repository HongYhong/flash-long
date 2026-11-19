"""Long-read exon & fuzzy splice-site corrector using reference gene profiles."""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

from isoquant_next.core.utils import contains_well_inside, equal_ranges, junctions_from_blocks

logger = logging.getLogger("IsoQuant")


class ExonCorrector:
    """Corrects fuzzy junctions, intron shifts, skipped exons, and fake terminal exons."""

    def __init__(self, gene_info, params, chr_record=None, profile_constructor=None) -> None:
        self.gene_info = gene_info
        self.chr_record = chr_record
        self.params = params
        self.delta = params.delta
        if profile_constructor is not None:
            self.intron_profile_constructor = profile_constructor
        else:
            from isoquant_next.assignment.profiles import OverlappingFeaturesProfileConstructor
            from functools import partial

            self.intron_profile_constructor = OverlappingFeaturesProfileConstructor(
                self.gene_info.intron_profiles.features,
                (self.gene_info.start, self.gene_info.end),
                comparator=partial(equal_ranges, delta=self.params.delta),
            )

    def correct_assigned_read(self, alignment_info, read_assignment) -> List[Tuple[int, int]]:
        if (
            len(alignment_info.read_exons) == 1
            or read_assignment.assignment_type.is_noninformative()
            or not read_assignment.isoform_matches
        ):
            return alignment_info.read_exons

        read_region, new_introns = self.correct_misalignments(alignment_info, read_assignment)
        if new_introns:
            corrected_exons = [(read_region[0], new_introns[0][0] - 1)]
            corrected_exons += junctions_from_blocks(new_introns)
            corrected_exons += [(new_introns[-1][1] + 1, read_region[1])]
        else:
            corrected_exons = [(read_region[0], read_region[1])]
        return corrected_exons

    def correct_fuzzy_junctions(self, alignment_info, read_assignment):
        intron_profile = alignment_info.combined_profile.read_intron_profile
        isoform_id = read_assignment.isoform_matches[0].assigned_transcript
        isoform_profile = self.gene_info.intron_profiles.profiles[isoform_id]
        ref_introns = self.gene_info.intron_profiles.features

        new_introns: List[Tuple[int, int]] = []
        read_pos, ref_pos = self.find_next(intron_profile, -1, -1)
        while read_pos < len(intron_profile.read_profile) and ref_pos < len(intron_profile.gene_profile):
            if isoform_profile[ref_pos] != 1:
                ref_pos += 1
            assert isoform_profile[ref_pos] == 1
            read_intron = intron_profile.read_features[read_pos]
            ref_intron = ref_introns[ref_pos]
            assert equal_ranges(ref_intron, read_intron, self.delta)
            new_introns.append(ref_intron)
            read_pos, ref_pos = self.find_next(intron_profile, read_pos, ref_pos)
        return (alignment_info.read_start, alignment_info.read_end), new_introns

    @staticmethod
    def find_next(intron_profile, read_pos: int, ref_pos: int) -> Tuple[int, int]:
        try:
            read_pos = intron_profile.read_profile.index(1, read_pos + 1)
            ref_pos = intron_profile.gene_profile.index(1, ref_pos + 1)
        except ValueError:
            read_pos = len(intron_profile.read_profile)
            ref_pos = len(intron_profile.gene_profile)
        return read_pos, ref_pos

    def correct_misalignments(self, alignment_info, read_assignment):
        from isoquant_next.assignment.types import MatchEventSubtype, SupplementaryMatchConstants

        event_map = {}
        for e in read_assignment.isoform_matches[0].match_subclassifications:
            if e.read_region == SupplementaryMatchConstants.undefined_region:
                continue
            if e.read_region[0] == SupplementaryMatchConstants.absent_position:
                if e.event_type == MatchEventSubtype.fake_micro_intron_retention and self.params.correct_microintron_retention:
                    event_map[-e.read_region[1] - 1] = e
            else:
                event_map[e.read_region[0]] = e

        intron_profile = alignment_info.combined_profile.read_intron_profile
        read_introns = intron_profile.read_features
        read_region = (alignment_info.read_start, alignment_info.read_end)
        isoform_id = read_assignment.isoform_matches[0].assigned_transcript
        isoform_region = self.gene_info.transcript_region(isoform_id)
        isoform_introns = self.gene_info.all_isoforms_introns[isoform_id]
        return self.process_events(alignment_info, event_map, read_region, read_introns, isoform_region, isoform_introns)

    def process_events(self, alignment_info, event_map, read_region, read_introns, isoform_region, isoform_introns):
        from isoquant_next.assignment.types import MatchEventSubtype

        if self.params.correct_fuzzy_junctions:
            potential_introns = self.intron_profile_constructor.match_genomic_features(read_introns)
            corrected_introns: List[Tuple[int, int]] = []
            assert len(potential_introns) == len(read_introns)
            for i in range(len(read_introns)):
                read_intron = read_introns[i]
                ref_intron = potential_introns[i]
                if read_intron[0] == ref_intron[0]:
                    left_site = read_intron[0]
                else:
                    start = min(read_intron[0], ref_intron[0])
                    end = max(read_intron[0], ref_intron[0]) - 1
                    indel_cnt, mm_cnt = alignment_info.get_error_count(
                        start, end, intron_index=i, left_site=True, chr_record=self.chr_record
                    )
                    left_site = read_intron[0] if indel_cnt == 0 and mm_cnt <= 1 else ref_intron[0]

                if read_intron[1] == ref_intron[1]:
                    right_site = read_intron[1]
                else:
                    start = min(read_intron[1], ref_intron[1]) + 1
                    end = max(read_intron[1], ref_intron[1])
                    indel_cnt, mm_cnt = alignment_info.get_error_count(
                        start, end, intron_index=i, left_site=False, chr_record=self.chr_record
                    )
                    right_site = read_intron[1] if indel_cnt == 0 and mm_cnt <= 1 else ref_intron[1]

                corrected_introns.append((left_site, right_site))
        else:
            corrected_introns = list(read_introns)

        corrected_read_region = read_region
        new_introns: List[Tuple[int, int]] = []
        i = 0
        misalignment_set = set()
        if self.params.correct_intron_shifts:
            misalignment_set.add(MatchEventSubtype.intron_shift)
        if self.params.correct_skipped_exons:
            misalignment_set.add(MatchEventSubtype.exon_misalignment)

        passthrough_events = {
            MatchEventSubtype.extra_intron_known,
            MatchEventSubtype.intron_alternation_known,
            MatchEventSubtype.intron_migration,
            MatchEventSubtype.exon_skipping_known,
            MatchEventSubtype.exon_merge_known,
            MatchEventSubtype.terminal_exon_shift_known,
            MatchEventSubtype.mutually_exclusive_exons_known,
            MatchEventSubtype.exon_gain_known,
            MatchEventSubtype.exon_detach_known,
            MatchEventSubtype.alternative_structure_known,
            MatchEventSubtype.alternative_structure_novel,
        }

        while i < len(corrected_introns):
            if -i - 1 in event_map:
                new_introns.append(isoform_introns[event_map[-i - 1].isoform_region[0]])
            if i not in event_map:
                new_introns.append(corrected_introns[i])
                i += 1
                continue

            event = event_map[i]
            if event.event_type == MatchEventSubtype.fake_terminal_exon_left and self.params.correct_fake_terminal_exons:
                corrected_read_region = (read_introns[event.read_region[0]][1] + 1, corrected_read_region[1])
            elif event.event_type == MatchEventSubtype.fake_terminal_exon_right and self.params.correct_fake_terminal_exons:
                corrected_read_region = (corrected_read_region[0], read_introns[event.read_region[0]][0] - 1)
            elif event.event_type == MatchEventSubtype.terminal_exon_misalignment_left and self.params.correct_terminal_exons:
                new_introns.append(isoform_introns[event.isoform_region[0]])
                corrected_read_region = (isoform_region[0], corrected_read_region[1])
            elif event.event_type == MatchEventSubtype.terminal_exon_misalignment_right and self.params.correct_terminal_exons:
                new_introns.append(isoform_introns[event.isoform_region[0]])
                corrected_read_region = (corrected_read_region[0], isoform_region[1])
            elif (
                event.event_type in misalignment_set
                and contains_well_inside(
                    read_region,
                    (isoform_introns[event.isoform_region[0]][0], isoform_introns[event.isoform_region[1]][1]),
                    self.params.delta,
                )
            ):
                new_introns.extend(isoform_introns[k] for k in range(event.isoform_region[0], event.isoform_region[1] + 1))
            elif event.event_type in passthrough_events:
                new_introns.extend(corrected_introns[k] for k in range(event.read_region[0], event.read_region[1] + 1))
            else:
                new_introns.extend(read_introns[k] for k in range(event.read_region[0], event.read_region[1] + 1))

            i = event.read_region[1] + 1

        return corrected_read_region, new_introns
