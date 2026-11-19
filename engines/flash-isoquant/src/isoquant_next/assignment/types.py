"""Read assignment enums, events, IsoformMatch, ReadAssignment, CompactAssignmentRecord, and PolyA verifier."""

from __future__ import annotations

import logging
import math
import struct
from enum import Enum, unique
from typing import Dict, List, Optional, Sequence, Tuple

from isoquant_next.alignment.alignment_info import PolyAInfo, shift_polya, shift_polyt
from isoquant_next.core.id_policy import SimpleIDDistributor
from isoquant_next.core.serialization import (
    SHORT_FLOAT_MULTIPLIER,
    SHORT_INT_BYTES,
    read_bool_array,
    read_dict,
    read_int,
    read_int_neg,
    read_int_or_none,
    read_list,
    read_list_of_pairs,
    read_short_int,
    read_string,
    write_bool_array,
    write_dict,
    write_int,
    write_int_neg,
    write_int_or_none,
    write_list,
    write_list_of_pairs,
    write_short_int,
    write_string,
)
from isoquant_next.core.utils import interval_len, intervals_total_length, junctions_from_blocks

logger = logging.getLogger("IsoQuant")


@unique
class ReadAssignmentType(Enum):
    unique = 1
    noninformative = 0
    intergenic = 20
    ambiguous = 10
    unique_minor_difference = 2
    inconsistent = 30
    inconsistent_non_intronic = 31
    inconsistent_ambiguous = 32
    inconsistent_genic = 33
    inconsistent_multigenic = 34
    discarded = 40
    suspended = 255

    def is_unique(self) -> bool:
        return self in (ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference)

    def is_ambiguous(self) -> bool:
        return self in (
            ReadAssignmentType.ambiguous,
            ReadAssignmentType.inconsistent_ambiguous,
            ReadAssignmentType.inconsistent_multigenic,
        )

    def is_consistent(self) -> bool:
        return self in (
            ReadAssignmentType.unique,
            ReadAssignmentType.unique_minor_difference,
            ReadAssignmentType.ambiguous,
        )

    def is_inconsistent(self) -> bool:
        return self in (
            ReadAssignmentType.inconsistent,
            ReadAssignmentType.inconsistent_non_intronic,
            ReadAssignmentType.inconsistent_ambiguous,
            ReadAssignmentType.inconsistent_genic,
            ReadAssignmentType.inconsistent_multigenic,
        )

    def is_noninformative(self) -> bool:
        return self in (ReadAssignmentType.noninformative, ReadAssignmentType.intergenic)

    def is_unassigned(self) -> bool:
        return self in (
            ReadAssignmentType.noninformative,
            ReadAssignmentType.intergenic,
            ReadAssignmentType.discarded,
            ReadAssignmentType.suspended,
        )

    def __lt__(self, other: "ReadAssignmentType") -> bool:
        return self.value < other.value


@unique
class MatchClassification(Enum):
    undefined = 0
    full_splice_match = 10
    incomplete_splice_match = 11
    mono_exon_match = 12
    novel_in_catalog = 20
    novel_not_in_catalog = 21
    genic = 30
    antisense = 40
    fusion = 50
    intergenic = 32
    genic_intron = 31

    @staticmethod
    def get_mono_exon_classification(match_event_subtypes) -> "MatchClassification":
        if any(
            e.event_type
            in {
                MatchEventSubtype.alternative_polya_site_left,
                MatchEventSubtype.alternative_polya_site_right,
                MatchEventSubtype.internal_polya_left,
                MatchEventSubtype.internal_polya_right,
            }
            for e in match_event_subtypes
        ):
            return MatchClassification.novel_not_in_catalog
        if any(e.event_type == MatchEventSubtype.unspliced_intron_retention for e in match_event_subtypes):
            return MatchClassification.novel_in_catalog
        if any(
            e.event_type
            in {
                MatchEventSubtype.incomplete_intron_retention_left,
                MatchEventSubtype.incomplete_intron_retention_right,
            }
            for e in match_event_subtypes
        ):
            return MatchClassification.genic
        if match_event_subtypes[0].event_type == MatchEventSubtype.fake_micro_intron_retention:
            return MatchClassification.incomplete_splice_match
        if match_event_subtypes[0].event_type == MatchEventSubtype.mono_exon_match:
            return MatchClassification.mono_exon_match
        if match_event_subtypes[0].event_type == MatchEventSubtype.mono_exonic:
            return MatchClassification.incomplete_splice_match
        logger.warning("Unexpected set of monoexonic events: %s", match_event_subtypes)
        return MatchClassification.undefined

    @staticmethod
    def get_inconsistency_classification(match_event_subtypes) -> "MatchClassification":
        if any(me.event_type in nnic_event_types for me in match_event_subtypes):
            return MatchClassification.novel_not_in_catalog
        if any(me.event_type in nic_event_types for me in match_event_subtypes):
            return MatchClassification.novel_in_catalog
        return MatchClassification.undefined


@unique
class MatchEventSubtype(Enum):
    none = 0
    undefined = 1
    mono_exonic = 10
    antisense = 9
    fsm = 11
    ism_left = 15
    ism_right = 13
    ism_internal = 14
    mono_exon_match = 19
    intron_shift = 21
    exon_misalignment = 22
    fake_terminal_exon_left = 24
    fake_terminal_exon_right = 26
    terminal_exon_misalignment_right = 28
    terminal_exon_misalignment_left = 29
    exon_elongation_left = 25
    exon_elongation_right = 23
    intron_retention = 31
    unspliced_intron_retention = 32
    incomplete_intron_retention_left = 38
    incomplete_intron_retention_right = 39
    fake_micro_intron_retention = 34
    alt_left_site_known = 101
    alt_right_site_known = 102
    alt_left_site_novel = 103
    alt_right_site_novel = 104
    extra_intron_novel = 1012
    extra_intron_known = 1011
    extra_intron_flanking_left = 1013
    extra_intron_flanking_right = 1015
    major_exon_elongation_left = 1005
    major_exon_elongation_right = 1003
    intron_migration = 114
    intron_alternation_novel = 115
    intron_alternation_known = 116
    mutually_exclusive_exons_novel = 121
    mutually_exclusive_exons_known = 122
    exon_skipping_known = 133
    exon_skipping_novel = 134
    exon_merge_known = 135
    exon_merge_novel = 136
    exon_gain_known = 145
    exon_gain_novel = 146
    exon_detach_known = 147
    exon_detach_novel = 148
    terminal_exon_shift_known = 151
    terminal_exon_shift_novel = 152
    alternative_structure_novel = 181
    alternative_structure_known = 182
    alternative_polya_site_right = 200
    alternative_polya_site_left = 201
    internal_polya_left = 202
    internal_polya_right = 203
    alternative_tss_left = 204
    alternative_tss_right = 205
    correct_polya_site_left = 222
    correct_polya_site_right = 223
    aligned_polya_tail = 300
    terminal_site_match_left = 250
    terminal_site_match_right = 251
    terminal_site_match_left_precise = 252
    terminal_site_match_right_precise = 253

    def __lt__(self, other: "MatchEventSubtype") -> bool:
        return self.value < other.value

    @staticmethod
    def is_alignment_artifact(match_event_subtype) -> bool:
        return match_event_subtype in {
            MatchEventSubtype.intron_shift,
            MatchEventSubtype.exon_misalignment,
            MatchEventSubtype.terminal_exon_misalignment_right,
            MatchEventSubtype.terminal_exon_misalignment_left,
            MatchEventSubtype.fake_terminal_exon_left,
            MatchEventSubtype.fake_terminal_exon_right,
            MatchEventSubtype.fake_micro_intron_retention,
        }

    @staticmethod
    def is_minor_error(match_event_subtype) -> bool:
        return match_event_subtype in {
            MatchEventSubtype.exon_elongation_left,
            MatchEventSubtype.exon_elongation_right,
            MatchEventSubtype.intron_shift,
            MatchEventSubtype.exon_misalignment,
            MatchEventSubtype.fake_micro_intron_retention,
            MatchEventSubtype.fake_terminal_exon_left,
            MatchEventSubtype.fake_terminal_exon_right,
            MatchEventSubtype.terminal_exon_misalignment_left,
            MatchEventSubtype.terminal_exon_misalignment_right,
        }

    @staticmethod
    def is_consistent(match_event_subtype) -> bool:
        return match_event_subtype in {
            MatchEventSubtype.none,
            MatchEventSubtype.mono_exonic,
            MatchEventSubtype.mono_exon_match,
            MatchEventSubtype.fsm,
            MatchEventSubtype.ism_left,
            MatchEventSubtype.ism_right,
            MatchEventSubtype.ism_internal,
            MatchEventSubtype.correct_polya_site_left,
            MatchEventSubtype.correct_polya_site_right,
            MatchEventSubtype.terminal_site_match_left,
            MatchEventSubtype.terminal_site_match_right,
            MatchEventSubtype.terminal_site_match_left_precise,
            MatchEventSubtype.terminal_site_match_right_precise,
        }

    @staticmethod
    def is_major_elongation(match_event_subtype) -> bool:
        return match_event_subtype in {
            MatchEventSubtype.major_exon_elongation_left,
            MatchEventSubtype.major_exon_elongation_right,
        }

    @staticmethod
    def is_minor_elongation(match_event_subtype) -> bool:
        return match_event_subtype in {
            MatchEventSubtype.exon_elongation_left,
            MatchEventSubtype.exon_elongation_right,
        }

    @staticmethod
    def is_major_inconsistency(match_event_subtype) -> bool:
        return match_event_subtype in all_major_events

    @staticmethod
    def is_intronic_inconsistency(match_event_subtype) -> bool:
        return match_event_subtype in intronic_major_events


event_subtype_cost: Dict[MatchEventSubtype, float] = {
    MatchEventSubtype.none: 0,
    MatchEventSubtype.undefined: 0,
    MatchEventSubtype.mono_exonic: 0,
    MatchEventSubtype.fsm: 0,
    MatchEventSubtype.ism_left: 0,
    MatchEventSubtype.ism_right: 0,
    MatchEventSubtype.ism_internal: 0,
    MatchEventSubtype.mono_exon_match: 0,
    MatchEventSubtype.intron_shift: 0.1,
    MatchEventSubtype.exon_misalignment: 0.1,
    MatchEventSubtype.terminal_exon_misalignment_right: 0.1,
    MatchEventSubtype.terminal_exon_misalignment_left: 0.1,
    MatchEventSubtype.fake_terminal_exon_left: 0.2,
    MatchEventSubtype.fake_terminal_exon_right: 0.2,
    MatchEventSubtype.exon_elongation_left: 0.1,
    MatchEventSubtype.exon_elongation_right: 0.1,
    MatchEventSubtype.intron_retention: 0.6,
    MatchEventSubtype.unspliced_intron_retention: 0.6,
    MatchEventSubtype.incomplete_intron_retention_left: 0.7,
    MatchEventSubtype.incomplete_intron_retention_right: 0.7,
    MatchEventSubtype.fake_micro_intron_retention: 0.1,
    MatchEventSubtype.alt_left_site_known: 1,
    MatchEventSubtype.alt_right_site_known: 1,
    MatchEventSubtype.alt_left_site_novel: 1,
    MatchEventSubtype.alt_right_site_novel: 1,
    MatchEventSubtype.extra_intron_novel: 1,
    MatchEventSubtype.extra_intron_known: 1,
    MatchEventSubtype.extra_intron_flanking_left: 1,
    MatchEventSubtype.extra_intron_flanking_right: 1,
    MatchEventSubtype.major_exon_elongation_left: 0.6,
    MatchEventSubtype.major_exon_elongation_right: 0.6,
    MatchEventSubtype.intron_migration: 1,
    MatchEventSubtype.intron_alternation_novel: 1,
    MatchEventSubtype.intron_alternation_known: 1,
    MatchEventSubtype.mutually_exclusive_exons_novel: 0.8,
    MatchEventSubtype.mutually_exclusive_exons_known: 0.8,
    MatchEventSubtype.exon_skipping_known: 1,
    MatchEventSubtype.exon_skipping_novel: 1,
    MatchEventSubtype.exon_merge_known: 0.75,
    MatchEventSubtype.exon_merge_novel: 0.75,
    MatchEventSubtype.exon_gain_known: 1,
    MatchEventSubtype.exon_gain_novel: 1,
    MatchEventSubtype.exon_detach_known: 0.75,
    MatchEventSubtype.exon_detach_novel: 0.75,
    MatchEventSubtype.terminal_exon_shift_known: 0.5,
    MatchEventSubtype.terminal_exon_shift_novel: 0.5,
    MatchEventSubtype.alternative_structure_novel: 1,
    MatchEventSubtype.alternative_structure_known: 1,
    MatchEventSubtype.alternative_polya_site_left: 0.6,
    MatchEventSubtype.alternative_polya_site_right: 0.6,
    MatchEventSubtype.internal_polya_left: 0.5,
    MatchEventSubtype.internal_polya_right: 0.5,
    MatchEventSubtype.alternative_tss_left: 0.6,
    MatchEventSubtype.alternative_tss_right: 0.6,
    MatchEventSubtype.correct_polya_site_left: 0,
    MatchEventSubtype.correct_polya_site_right: 0,
    MatchEventSubtype.aligned_polya_tail: 0,
    MatchEventSubtype.terminal_site_match_left: 0,
    MatchEventSubtype.terminal_site_match_right: 0,
    MatchEventSubtype.terminal_site_match_left_precise: 0,
    MatchEventSubtype.terminal_site_match_right_precise: 0,
}


def elongation_cost(params, elongation_len: int) -> float:
    min_cost = event_subtype_cost[MatchEventSubtype.exon_elongation_left]
    max_cost = event_subtype_cost[MatchEventSubtype.major_exon_elongation_left]
    lower_bound = params.minor_exon_extension
    upper_bound = params.major_exon_extension
    if elongation_len <= lower_bound:
        return min_cost
    if elongation_len >= upper_bound:
        return max_cost
    return min_cost + (max_cost - min_cost) * (elongation_len - lower_bound) / (upper_bound - lower_bound)


nnic_event_types = {
    MatchEventSubtype.alt_left_site_novel,
    MatchEventSubtype.alt_right_site_novel,
    MatchEventSubtype.extra_intron_novel,
    MatchEventSubtype.extra_intron_flanking_left,
    MatchEventSubtype.extra_intron_flanking_right,
    MatchEventSubtype.mutually_exclusive_exons_novel,
    MatchEventSubtype.exon_gain_novel,
    MatchEventSubtype.exon_skipping_novel,
    MatchEventSubtype.exon_detach_novel,
    MatchEventSubtype.exon_merge_novel,
    MatchEventSubtype.terminal_exon_shift_novel,
    MatchEventSubtype.alternative_structure_novel,
    MatchEventSubtype.intron_alternation_novel,
}

nic_event_types = {
    MatchEventSubtype.unspliced_intron_retention,
    MatchEventSubtype.intron_retention,
    MatchEventSubtype.alt_left_site_known,
    MatchEventSubtype.alt_right_site_known,
    MatchEventSubtype.extra_intron_known,
    MatchEventSubtype.intron_migration,
    MatchEventSubtype.mutually_exclusive_exons_known,
    MatchEventSubtype.exon_skipping_known,
    MatchEventSubtype.exon_detach_known,
    MatchEventSubtype.exon_merge_known,
    MatchEventSubtype.terminal_exon_shift_known,
    MatchEventSubtype.exon_gain_known,
    MatchEventSubtype.alternative_structure_known,
    MatchEventSubtype.intron_alternation_known,
    MatchEventSubtype.major_exon_elongation_left,
    MatchEventSubtype.major_exon_elongation_right,
    MatchEventSubtype.incomplete_intron_retention_left,
    MatchEventSubtype.incomplete_intron_retention_right,
    MatchEventSubtype.internal_polya_right,
    MatchEventSubtype.internal_polya_left,
    MatchEventSubtype.alternative_polya_site_left,
    MatchEventSubtype.alternative_polya_site_right,
    MatchEventSubtype.alternative_tss_right,
    MatchEventSubtype.alternative_tss_left,
}

nonintronic_events = {
    MatchEventSubtype.alternative_polya_site_left,
    MatchEventSubtype.alternative_polya_site_right,
    MatchEventSubtype.alternative_tss_left,
    MatchEventSubtype.alternative_tss_right,
    MatchEventSubtype.internal_polya_right,
    MatchEventSubtype.internal_polya_left,
    MatchEventSubtype.major_exon_elongation_left,
    MatchEventSubtype.major_exon_elongation_right,
    MatchEventSubtype.exon_elongation_left,
    MatchEventSubtype.exon_elongation_right,
}

all_major_events = nic_event_types.union(nnic_event_types)
intronic_major_events = all_major_events.difference(nonintronic_events)

alternative_sites = {
    ("left", True): MatchEventSubtype.alt_left_site_known,
    ("left", False): MatchEventSubtype.alt_left_site_novel,
    ("right", True): MatchEventSubtype.alt_right_site_known,
    ("right", False): MatchEventSubtype.alt_right_site_novel,
}


class SupplementaryMatchConstants:
    extra_left_mod_position = (1 << 30) - 1
    extra_right_mod_position = (1 << 30) + 1
    undefined_position = 1 << 31
    undefined_region = (undefined_position, undefined_position)
    extra_left_region = (extra_left_mod_position, extra_left_mod_position)
    extra_right_region = (extra_right_mod_position, extra_right_mod_position)
    absent_position = (1 << 31) - 1


class MatchEvent:
    __slots__ = ("event_type", "isoform_region", "read_region", "event_info")

    def __init__(
        self,
        event_type: MatchEventSubtype,
        isoform_region: tuple = SupplementaryMatchConstants.undefined_region,
        read_region: tuple = SupplementaryMatchConstants.undefined_region,
        event_info: int = 0,
    ) -> None:
        self.event_type = event_type
        self.isoform_region = isoform_region
        self.read_region = read_region
        self.event_info = event_info

    @classmethod
    def deserialize(cls, infile) -> "MatchEvent":
        ev = cls.__new__(cls)
        ev.event_type = MatchEventSubtype(read_int(infile, SHORT_INT_BYTES))
        ev.isoform_region = (read_int(infile), read_int(infile))
        ev.read_region = (read_int(infile), read_int(infile))
        ev.event_info = read_int_neg(infile)
        return ev

    def serialize(self, outfile) -> None:
        write_int(self.event_type.value, outfile, SHORT_INT_BYTES)
        write_int(self.isoform_region[0], outfile)
        write_int(self.isoform_region[1], outfile)
        write_int(self.read_region[0], outfile)
        write_int(self.read_region[1], outfile)
        write_int_neg(self.event_info, outfile)

    def __repr__(self) -> str:
        return f"{self.event_type.name}:{self.isoform_region},{self.read_region},{self.event_info}"


class IsoformMatch:
    """Memory-compact isoform match storing interned gene/transcript integer IDs."""

    __slots__ = (
        "_string_pools",
        "assigned_gene_id",
        "assigned_transcript_id",
        "transcript_strand",
        "match_classification",
        "penalty_score",
        "match_subclassifications",
    )

    def __init__(
        self,
        match_classification: MatchClassification,
        string_pools,
        assigned_gene: Optional[str] = None,
        assigned_transcript: Optional[str] = None,
        match_subclassification=None,
        transcript_strand: str = ".",
        penalty_score: float = 0.0,
    ) -> None:
        assert string_pools is not None, "string_pools is required"
        self._string_pools = string_pools
        self.assigned_gene_id: Optional[int] = string_pools.gene_pool.get_int(assigned_gene) if assigned_gene else None
        self.assigned_transcript_id: Optional[int] = (
            string_pools.transcript_pool.get_int(assigned_transcript) if assigned_transcript else None
        )
        self.transcript_strand = transcript_strand
        self.match_classification = match_classification
        self.penalty_score = penalty_score
        if match_subclassification is None:
            self.match_subclassifications: List[MatchEvent] = []
        elif isinstance(match_subclassification, list):
            self.match_subclassifications = [x for x in match_subclassification if x.event_type != MatchEventSubtype.none]
        else:
            self.match_subclassifications = [match_subclassification]

    @property
    def assigned_gene(self) -> Optional[str]:
        return self._string_pools.gene_pool.get_str(self.assigned_gene_id) if self.assigned_gene_id is not None else None

    @assigned_gene.setter
    def assigned_gene(self, value: Optional[str]) -> None:
        self.assigned_gene_id = self._string_pools.gene_pool.get_int(value) if value is not None else None

    @property
    def assigned_transcript(self) -> Optional[str]:
        return (
            self._string_pools.transcript_pool.get_str(self.assigned_transcript_id)
            if self.assigned_transcript_id is not None
            else None
        )

    @assigned_transcript.setter
    def assigned_transcript(self, value: Optional[str]) -> None:
        self.assigned_transcript_id = self._string_pools.transcript_pool.get_int(value) if value is not None else None

    @classmethod
    def deserialize(cls, infile, string_pools) -> "IsoformMatch":
        assert string_pools is not None, "string_pools is required"
        m = cls.__new__(cls)
        m._string_pools = string_pools
        m.assigned_gene_id = read_int_or_none(infile)
        m.assigned_transcript_id = read_int_or_none(infile)
        m.transcript_strand = read_string(infile)
        m.match_classification = MatchClassification(read_short_int(infile))
        m.penalty_score = float(read_int(infile)) / float(SHORT_FLOAT_MULTIPLIER)
        m.match_subclassifications = read_list(infile, MatchEvent.deserialize)
        return m

    def serialize(self, outfile) -> None:
        write_int_or_none(self.assigned_gene_id, outfile)
        write_int_or_none(self.assigned_transcript_id, outfile)
        write_string(self.transcript_strand, outfile)
        write_short_int(self.match_classification.value, outfile)
        write_int(int(self.penalty_score * SHORT_FLOAT_MULTIPLIER), outfile)
        write_list(self.match_subclassifications, outfile, MatchEvent.serialize)

    def add_subclassification(self, match_subclassification: MatchEvent) -> None:
        if len(self.match_subclassifications) == 1 and self.match_subclassifications[0].event_type in {
            MatchEventSubtype.undefined,
            MatchEventSubtype.none,
        }:
            self.match_subclassifications = [match_subclassification]
        else:
            self.match_subclassifications.append(match_subclassification)

    def set_classification(self, classification: MatchClassification) -> None:
        self.match_classification = classification

    def monoexon_is_consistent(self) -> bool:
        valid = {MatchEventSubtype.none, MatchEventSubtype.mono_exonic, MatchEventSubtype.mono_exon_match}
        return all(el.event_type in valid for el in self.match_subclassifications)


# Fixed-width packed binary record replacing full ReadAssignment IPC pickling across workers
_COMPACT_HEADER_STRUCT = struct.Struct("<QiiiiBBHHd")


class CompactMultimapTuple:
    """Zero-pickling fixed-width record (binary packable via struct) for cross-process multimap resolution."""

    __slots__ = (
        "assignment_id",
        "read_id",
        "chr_id_int",
        "start",
        "end",
        "genomic_start",
        "genomic_end",
        "multimapper",
        "polyA_found",
        "assignment_type",
        "gene_assignment_type",
        "penalty_score",
        "gene_ids",
        "isoform_ids",
    )

    def __init__(
        self,
        assignment_id: int,
        read_id: str,
        chr_id_int: int,
        start: int,
        end: int,
        genomic_start: int,
        genomic_end: int,
        multimapper: bool,
        polyA_found: bool,
        assignment_type: ReadAssignmentType,
        gene_assignment_type: ReadAssignmentType,
        penalty_score: float,
        gene_ids: Tuple[int, ...],
        isoform_ids: Tuple[int, ...],
    ) -> None:
        self.assignment_id = assignment_id
        self.read_id = read_id
        self.chr_id_int = chr_id_int
        self.start = start
        self.end = end
        self.genomic_start = genomic_start
        self.genomic_end = genomic_end
        self.multimapper = multimapper
        self.polyA_found = polyA_found
        self.assignment_type = assignment_type
        self.gene_assignment_type = gene_assignment_type
        self.penalty_score = penalty_score
        self.gene_ids = gene_ids
        self.isoform_ids = isoform_ids

    @property
    def genomic_region(self) -> Tuple[int, int]:
        return (self.genomic_start, self.genomic_end)

    def to_bytes(self) -> bytes:
        rid_b = self.read_id.encode("utf-8")
        fixed = _COMPACT_HEADER_STRUCT.pack(
            self.assignment_id,
            self.chr_id_int,
            self.start,
            self.end,
            self.genomic_start,
            1 if self.multimapper else 0,
            1 if self.polyA_found else 0,
            self.assignment_type.value,
            self.gene_assignment_type.value,
            self.penalty_score,
        )
        n_g = len(self.gene_ids)
        n_t = len(self.isoform_ids)
        tail = struct.pack(f"<iHH{len(rid_b)}s{n_g}i{n_t}i", self.genomic_end, n_g, n_t, rid_b, *self.gene_ids, *self.isoform_ids)
        return fixed + tail

    @classmethod
    def from_bytes(cls, data: bytes) -> "CompactMultimapTuple":
        aid, cid, st, en, gs, mm, pa, atv, gtv, pen = _COMPACT_HEADER_STRUCT.unpack_from(data, 0)
        off = _COMPACT_HEADER_STRUCT.size
        ge, n_g, n_t = struct.unpack_from("<iHH", data, off)
        off += 8
        # Compute rid length from remaining payload bytes
        rid_len = len(data) - off - 4 * (n_g + n_t)
        rid = data[off : off + rid_len].decode("utf-8")
        off += rid_len
        g_ids = struct.unpack_from(f"<{n_g}i", data, off) if n_g else ()
        off += 4 * n_g
        t_ids = struct.unpack_from(f"<{n_t}i", data, off) if n_t else ()
        return cls(
            aid,
            rid,
            cid,
            st,
            en,
            gs,
            ge,
            bool(mm),
            bool(pa),
            ReadAssignmentType(atv),
            ReadAssignmentType(gtv),
            pen,
            g_ids,
            t_ids,
        )


class BasicReadAssignment:
    """Compact slots representation for serialization and multimap resolution."""

    __slots__ = (
        "assignment_id",
        "read_id",
        "_string_pools",
        "chr_id_int",
        "barcode_id",
        "umi_id",
        "start",
        "end",
        "genomic_region",
        "multimapper",
        "polyA_found",
        "assignment_type",
        "gene_assignment_type",
        "penalty_score",
        "gene_ids",
        "isoform_ids",
        "read_group_ids",
        "_chr_id_str",
        "_genes",
        "_isoforms",
    )

    def __init__(self, read_assignment: "ReadAssignment") -> None:
        self.assignment_id = read_assignment.assignment_id
        self.read_id = read_assignment.read_id
        assert hasattr(read_assignment, "_string_pools") and read_assignment._string_pools is not None
        self._string_pools = read_assignment._string_pools
        self.chr_id_int = read_assignment.chr_id_int
        self.barcode_id = read_assignment.barcode_id
        self.umi_id = read_assignment.umi_id
        self.read_group_ids = list(read_assignment.read_group_ids)
        self.start = 0
        self.end = 0
        if read_assignment.exons:
            self.start = read_assignment.exons[0][0]
            self.end = read_assignment.exons[-1][1]
        self.genomic_region = read_assignment.genomic_region
        self.multimapper = read_assignment.multimapper
        self.polyA_found = read_assignment.polyA_found
        self.assignment_type = read_assignment.assignment_type
        self.gene_assignment_type = read_assignment.gene_assignment_type
        self.penalty_score = 0.0
        self.gene_ids: List[int] = []
        self.isoform_ids: List[int] = []
        if read_assignment.isoform_matches:
            gene_set = set()
            iso_set = set()
            for m in read_assignment.isoform_matches:
                self.penalty_score = min(self.penalty_score, read_assignment.isoform_matches[0].penalty_score)
                if m.assigned_gene_id is not None:
                    gene_set.add(m.assigned_gene_id)
                if m.assigned_transcript_id is not None:
                    iso_set.add(m.assigned_transcript_id)
            self.gene_ids = list(gene_set)
            self.isoform_ids = list(iso_set)

    def to_compact_tuple(self) -> CompactMultimapTuple:
        return CompactMultimapTuple(
            self.assignment_id,
            self.read_id,
            self.chr_id_int if self.chr_id_int is not None else -1,
            self.start,
            self.end,
            self.genomic_region[0],
            self.genomic_region[1],
            self.multimapper,
            self.polyA_found,
            self.assignment_type,
            self.gene_assignment_type,
            self.penalty_score,
            tuple(self.gene_ids),
            tuple(self.isoform_ids),
        )

    @property
    def chr_id(self) -> str:
        if self.chr_id_int is not None and self._string_pools is not None:
            return self._string_pools.chromosome_pool.get_str(self.chr_id_int)
        return getattr(self, "_chr_id_str", ".")

    @chr_id.setter
    def chr_id(self, value: Optional[str]) -> None:
        self.chr_id_int = self._string_pools.chromosome_pool.get_int(value) if value else None

    @property
    def barcode(self) -> Optional[str]:
        return self._string_pools.barcode_pool.get_str(self.barcode_id) if self.barcode_id is not None else None

    @barcode.setter
    def barcode(self, value: Optional[str]) -> None:
        self.barcode_id = self._string_pools.barcode_pool.get_int(value) if value else None

    @property
    def umi(self) -> Optional[str]:
        return self._string_pools.umi_pool.get_str(self.umi_id) if self.umi_id is not None else None

    @umi.setter
    def umi(self, value: Optional[str]) -> None:
        self.umi_id = self._string_pools.umi_pool.get_int(value) if value else None

    @property
    def genes(self) -> List[str]:
        if self.gene_ids and self._string_pools is not None:
            return [self._string_pools.gene_pool.get_str(gid) for gid in self.gene_ids]
        return getattr(self, "_genes", [])

    @genes.setter
    def genes(self, value: Sequence[str]) -> None:
        self.gene_ids = [self._string_pools.gene_pool.get_int(g) for g in value]

    @property
    def isoforms(self) -> List[str]:
        if self.isoform_ids and self._string_pools is not None:
            return [self._string_pools.transcript_pool.get_str(tid) for tid in self.isoform_ids]
        return getattr(self, "_isoforms", [])

    @isoforms.setter
    def isoforms(self, value: Sequence[str]) -> None:
        self.isoform_ids = [self._string_pools.transcript_pool.get_int(t) for t in value]

    def __eq__(self, other: object) -> bool:
        if isinstance(other, BasicReadAssignment):
            return (
                self.read_id == other.read_id
                and self.chr_id == other.chr_id
                and self.start == other.start
                and self.end == other.end
                and self.isoforms == other.isoforms
            )
        return False

    @classmethod
    def deserialize(cls, infile, string_pools) -> "BasicReadAssignment":
        assert string_pools is not None, "string_pools is required"
        ra = cls.__new__(cls)
        ra._string_pools = string_pools
        ra.assignment_id = read_int(infile)
        ra.read_id = read_string(infile)
        ra.chr_id_int = read_int_or_none(infile)
        ra.start = read_int(infile)
        ra.end = read_int(infile)
        ra.genomic_region = (read_int(infile), read_int(infile))
        bool_arr = read_bool_array(infile, 2)
        ra.multimapper = bool_arr[0]
        ra.polyA_found = bool_arr[1]
        ra.assignment_type = ReadAssignmentType(read_short_int(infile))
        ra.gene_assignment_type = ReadAssignmentType(read_short_int(infile))
        ra.penalty_score = float(read_int(infile)) / float(SHORT_FLOAT_MULTIPLIER)
        ra.gene_ids = read_list(infile, read_int)
        ra.isoform_ids = read_list(infile, read_int)
        ra.barcode_id = None
        ra.umi_id = None
        ra.read_group_ids = []
        return ra

    @classmethod
    def deserialize_from_read_assignment(cls, infile, string_pools) -> "BasicReadAssignment":
        assert string_pools is not None, "string_pools is required"
        ra = cls.__new__(cls)
        ra._string_pools = string_pools
        ra.assignment_id = read_int(infile)
        ra.read_id = read_string(infile)
        ra.genomic_region = (read_int(infile), read_int(infile))
        exons = read_list_of_pairs(infile, read_int)
        ra.start = exons[0][0]
        ra.end = exons[-1][1]
        read_list_of_pairs(infile, read_int)
        read_list_of_pairs(infile, read_int)
        bool_arr = read_bool_array(infile, 3)
        ra.multimapper = bool_arr[0]
        ra.polyA_found = bool_arr[1]
        read_int_neg(infile)
        read_int_neg(infile)
        read_int_neg(infile)
        read_int_neg(infile)
        ra.read_group_ids = read_list(infile, read_int)
        ra.barcode_id = read_int_or_none(infile)
        ra.umi_id = read_int_or_none(infile)
        read_string(infile)
        read_string(infile)
        ra.chr_id_int = read_int_or_none(infile)
        read_short_int(infile)
        ra.assignment_type = ReadAssignmentType(read_short_int(infile))
        ra.gene_assignment_type = ReadAssignmentType(read_short_int(infile))
        ra.penalty_score = 0.0
        isoform_matches = read_list(infile, lambda f: IsoformMatch.deserialize(f, string_pools))
        gene_set = set()
        iso_set = set()
        for m in isoform_matches:
            ra.penalty_score = min(ra.penalty_score, isoform_matches[0].penalty_score)
            if m.assigned_gene_id is not None:
                gene_set.add(m.assigned_gene_id)
            if m.assigned_transcript_id is not None:
                iso_set.add(m.assigned_transcript_id)
        ra.gene_ids = list(gene_set)
        ra.isoform_ids = list(iso_set)
        read_dict(infile)
        read_dict(infile)
        read_short_int(infile)
        read_list(infile, read_int_neg)
        read_list(infile, read_int_neg)
        return ra

    def serialize(self, outfile) -> None:
        write_int(self.assignment_id, outfile)
        write_string(self.read_id, outfile)
        write_int_or_none(self.chr_id_int, outfile)
        write_int(self.start, outfile)
        write_int(self.end, outfile)
        write_int(self.genomic_region[0], outfile)
        write_int(self.genomic_region[1], outfile)
        write_bool_array([self.multimapper, self.polyA_found], outfile)
        write_short_int(self.assignment_type.value, outfile)
        write_short_int(self.gene_assignment_type.value, outfile)
        write_int(int(self.penalty_score * SHORT_FLOAT_MULTIPLIER), outfile)
        write_list(self.gene_ids, outfile, write_int)
        write_list(self.isoform_ids, outfile, write_int)


class ReadAssignment:
    """Memory-compact ReadAssignment using __slots__ and StringPoolManager interned IDs."""

    assignment_id_generator = SimpleIDDistributor()

    __slots__ = (
        "assignment_id",
        "read_id",
        "genomic_region",
        "exons",
        "boundary_deletions",
        "corrected_exons",
        "corrected_introns",
        "gene_info",
        "multimapper",
        "polyA_found",
        "cage_found",
        "polya_info",
        "_string_pools",
        "read_group_ids",
        "barcode_id",
        "umi_id",
        "chr_id_int",
        "mapped_strand",
        "strand",
        "mapping_quality",
        "assignment_type",
        "isoform_matches",
        "gene_assignment_type",
        "additional_info",
        "additional_attributes",
        "introns_match",
        "exon_gene_profile",
        "intron_gene_profile",
    )

    def __init__(self, read_id: str, assignment_type: ReadAssignmentType, string_pools, match=None) -> None:
        assert string_pools is not None, "string_pools is required"
        self.assignment_id = ReadAssignment.assignment_id_generator.increment()
        self.read_id = read_id
        self.genomic_region = (0, 0)
        self.exons: Optional[List[Tuple[int, int]]] = None
        self.boundary_deletions = None
        self.corrected_exons: Optional[List[Tuple[int, int]]] = None
        self.corrected_introns: Optional[List[Tuple[int, int]]] = None
        self.gene_info = None
        self.multimapper = False
        self.polyA_found = False
        self.cage_found = False
        self.polya_info: Optional[PolyAInfo] = None

        self._string_pools = string_pools
        self.read_group_ids: List[int] = []
        self.barcode_id: Optional[int] = None
        self.umi_id: Optional[int] = None
        self.chr_id_int: Optional[int] = None

        self.mapped_strand = "."
        self.strand = "."
        self.mapping_quality = 0
        self.assignment_type = assignment_type
        if match is None:
            self.isoform_matches: List[IsoformMatch] = []
        elif isinstance(match, list):
            self.isoform_matches = match
        else:
            self.isoform_matches = [match]

        if self.assignment_type == ReadAssignmentType.ambiguous:
            assigned_genes = {m.assigned_gene for m in self.isoform_matches}
            self.gene_assignment_type = (
                ReadAssignmentType.ambiguous if len(assigned_genes) > 1 else ReadAssignmentType.unique
            )
        elif self.assignment_type == ReadAssignmentType.inconsistent_ambiguous:
            assigned_genes = {m.assigned_gene for m in self.isoform_matches}
            self.gene_assignment_type = (
                ReadAssignmentType.inconsistent_ambiguous
                if len(assigned_genes) > 1
                else ReadAssignmentType.inconsistent
            )
        elif self.assignment_type.is_unassigned():
            assigned_genes = {m.assigned_gene for m in self.isoform_matches if m.assigned_gene}
            if len(assigned_genes) > 1:
                self.gene_assignment_type = ReadAssignmentType.inconsistent_multigenic
            elif len(assigned_genes) == 1:
                self.gene_assignment_type = ReadAssignmentType.inconsistent_genic
            else:
                self.gene_assignment_type = self.assignment_type
        else:
            self.gene_assignment_type = self.assignment_type

        self.additional_info: Dict = {}
        self.additional_attributes: Dict = {}
        self.introns_match = False
        self.exon_gene_profile = []
        self.intron_gene_profile = []

    @property
    def chr_id(self) -> str:
        return self._string_pools.chromosome_pool.get_str(self.chr_id_int) if self.chr_id_int is not None else "."

    @chr_id.setter
    def chr_id(self, value: Optional[str]) -> None:
        self.chr_id_int = self._string_pools.chromosome_pool.get_int(value) if value else None

    @property
    def barcode(self) -> Optional[str]:
        return self._string_pools.barcode_pool.get_str(self.barcode_id) if self.barcode_id is not None else None

    @barcode.setter
    def barcode(self, value: Optional[str]) -> None:
        self.barcode_id = self._string_pools.barcode_pool.get_int(value) if value else None

    @property
    def umi(self) -> Optional[str]:
        return self._string_pools.umi_pool.get_str(self.umi_id) if self.umi_id is not None else None

    @umi.setter
    def umi(self, value: Optional[str]) -> None:
        self.umi_id = self._string_pools.umi_pool.get_int(value) if value else None

    @property
    def read_group(self) -> List[str]:
        return self._string_pools.read_group_from_ids(self.read_group_ids) if self.read_group_ids else []

    @read_group.setter
    def read_group(self, value: Optional[Sequence[str]]) -> None:
        self.read_group_ids = self._string_pools.read_group_to_ids(value) if value else []

    @classmethod
    def deserialize(cls, infile, gene_info, string_pools) -> "ReadAssignment":
        assert string_pools is not None, "string_pools is required"
        ra = cls.__new__(cls)
        ra._string_pools = string_pools
        ra.assignment_id = read_int(infile)
        ra.read_id = read_string(infile)
        ra.genomic_region = (read_int(infile), read_int(infile))
        ra.exons = read_list_of_pairs(infile, read_int)
        ra.boundary_deletions = read_list_of_pairs(infile, read_int)
        ra.corrected_exons = read_list_of_pairs(infile, read_int)
        ra.corrected_introns = junctions_from_blocks(ra.corrected_exons)
        ra.gene_info = gene_info
        bool_arr = read_bool_array(infile, 3)
        ra.multimapper = bool_arr[0]
        ra.polyA_found = bool_arr[1]
        ra.cage_found = bool_arr[2]
        ra.polya_info = PolyAInfo(
            read_int_neg(infile), read_int_neg(infile), read_int_neg(infile), read_int_neg(infile)
        )
        ra.read_group_ids = read_list(infile, read_int)
        ra.barcode_id = read_int_or_none(infile)
        ra.umi_id = read_int_or_none(infile)
        ra.mapped_strand = read_string(infile)
        ra.strand = read_string(infile)
        ra.chr_id_int = read_int_or_none(infile)
        ra.mapping_quality = read_short_int(infile)
        ra.assignment_type = ReadAssignmentType(read_short_int(infile))
        ra.gene_assignment_type = ReadAssignmentType(read_short_int(infile))
        ra.isoform_matches = read_list(infile, lambda f: IsoformMatch.deserialize(f, string_pools))
        ra.additional_info = read_dict(infile)
        ra.additional_attributes = read_dict(infile)
        ra.introns_match = bool(read_short_int(infile))
        ra.exon_gene_profile = read_list(infile, read_int_neg)
        ra.intron_gene_profile = read_list(infile, read_int_neg)
        return ra

    def serialize(self, outfile) -> None:
        pi = self.polya_info or PolyAInfo(-1, -1, -1, -1)
        write_int(self.assignment_id, outfile)
        write_string(self.read_id, outfile)
        write_int(self.genomic_region[0], outfile)
        write_int(self.genomic_region[1], outfile)
        write_list_of_pairs(self.exons or [], outfile, write_int)
        write_list_of_pairs(self.boundary_deletions or [], outfile, write_int)
        write_list_of_pairs(self.corrected_exons or [], outfile, write_int)
        write_bool_array([self.multimapper, self.polyA_found, self.cage_found], outfile)
        write_int_neg(pi.external_polya_pos, outfile)
        write_int_neg(pi.external_polyt_pos, outfile)
        write_int_neg(pi.internal_polya_pos, outfile)
        write_int_neg(pi.internal_polyt_pos, outfile)
        write_list(self.read_group_ids, outfile, write_int)
        write_int_or_none(self.barcode_id, outfile)
        write_int_or_none(self.umi_id, outfile)
        write_string(self.mapped_strand, outfile)
        write_string(self.strand, outfile)
        write_int_or_none(self.chr_id_int, outfile)
        write_short_int(self.mapping_quality, outfile)
        write_short_int(self.assignment_type.value, outfile)
        write_short_int(self.gene_assignment_type.value, outfile)
        write_list(self.isoform_matches, outfile, IsoformMatch.serialize)
        write_dict(self.additional_info, outfile)
        write_dict(self.additional_attributes, outfile)
        write_short_int(int(self.introns_match), outfile)
        write_list(self.exon_gene_profile, outfile, write_int_neg)
        write_list(self.intron_gene_profile, outfile, write_int_neg)

    def add_match(self, match: IsoformMatch) -> None:
        self.isoform_matches.append(match)

    def set_assignment_type(self, assignment_type: ReadAssignmentType) -> None:
        self.assignment_type = assignment_type

    def start(self) -> int:
        return self.exons[0][0]

    def end(self) -> int:
        return self.exons[-1][1]

    def length(self) -> int:
        return sum(x[1] - x[0] + 1 for x in (self.exons or ()))

    def exon_count(self) -> int:
        return len(self.exons or ())

    def set_additional_info(self, key, value) -> None:
        self.additional_info[key] = value

    def set_additional_attribute(self, key, value) -> None:
        self.additional_attributes[key] = value

    def add_match_attribute(self, match_event: MatchEvent) -> None:
        for m in self.isoform_matches:
            m.add_subclassification(match_event)


match_subtype_printable_names = {
    MatchEventSubtype.ism_left: ("ism_5", "ism_3", "ism"),
    MatchEventSubtype.ism_right: ("ism_3", "ism_5", "ism"),
    MatchEventSubtype.exon_elongation_left: ("exon_elongation_5", "exon_elongation_3", "exon_elongation"),
    MatchEventSubtype.exon_elongation_right: ("exon_elongation_3", "exon_elongation_5", "exon_elongation"),
    MatchEventSubtype.major_exon_elongation_left: ("major_exon_elongation_5", "major_exon_elongation_3", "major_exon_elongation"),
    MatchEventSubtype.major_exon_elongation_right: ("major_exon_elongation_3", "major_exon_elongation_5", "major_exon_elongation"),
    MatchEventSubtype.fake_terminal_exon_left: ("fake_terminal_exon_5", "fake_terminal_exon_3", "fake_terminal_exon"),
    MatchEventSubtype.fake_terminal_exon_right: ("fake_terminal_exon_3", "fake_terminal_exon_5", "fake_terminal_exon"),
    MatchEventSubtype.terminal_exon_misalignment_left: ("terminal_exon_misalignment_5", "terminal_exon_misalignment_3", "terminal_exon_misalignment"),
    MatchEventSubtype.terminal_exon_misalignment_right: ("terminal_exon_misalignment_3", "terminal_exon_misalignment_5", "terminal_exon_misalignment"),
    MatchEventSubtype.incomplete_intron_retention_left: ("incomplete_intron_retention_5", "incomplete_intron_retention_3", "incomplete_intron_retention"),
    MatchEventSubtype.incomplete_intron_retention_right: ("incomplete_intron_retention_3", "incomplete_intron_retention_5", "incomplete_intron_retention"),
    MatchEventSubtype.extra_intron_flanking_left: ("extra_intron_5", "extra_intron_3", "extra_intron"),
    MatchEventSubtype.extra_intron_flanking_right: ("extra_intron_3", "extra_intron_5", "extra_intron"),
    MatchEventSubtype.alt_left_site_known: ("alt_donor_site_known", "alt_acceptor_site_known", "alternative_splice_site_known"),
    MatchEventSubtype.alt_right_site_known: ("alt_acceptor_site_known", "alt_donor_site_known", "alternative_splice_site_known"),
    MatchEventSubtype.alt_left_site_novel: ("alt_donor_site_novel", "alt_acceptor_site_novel", "alternative_splice_site_novel"),
    MatchEventSubtype.alt_right_site_novel: ("alt_acceptor_site_novel", "alt_donor_site_novel", "alternative_splice_site_novel"),
    MatchEventSubtype.terminal_site_match_left: ("tss_match", "tes_match", "terminal_position_match"),
    MatchEventSubtype.terminal_site_match_right: ("tes_match", "tss_match", "terminal_position_match"),
    MatchEventSubtype.terminal_site_match_left_precise: ("tss_match_precise", "tes_match_precise", "terminal_position_match_precise"),
    MatchEventSubtype.terminal_site_match_right_precise: ("tes_match_precise", "tss_match_precise", "terminal_position_match_precise"),
}


def match_subtype_to_str(event: MatchEvent, strand: str) -> str:
    es = event.event_type
    if es in match_subtype_printable_names:
        if strand == "-":
            return match_subtype_printable_names[es][1]
        if strand == "+":
            return match_subtype_printable_names[es][0]
        return match_subtype_printable_names[es][2]
    return es.name


def regions_to_str(regions: Sequence[Tuple[int, int]]) -> str:
    return ",".join(f"{x[0]}-{x[1]}" for x in regions)


def match_subtype_to_str_with_additional_info(
    event: MatchEvent, strand: str, read_introns: Sequence[Tuple[int, int]], isoform_introns: Sequence[Tuple[int, int]]
) -> str:
    es = event.event_type
    add = ""
    if es in {
        MatchEventSubtype.intron_retention,
        MatchEventSubtype.unspliced_intron_retention,
        MatchEventSubtype.incomplete_intron_retention_left,
        MatchEventSubtype.incomplete_intron_retention_right,
        MatchEventSubtype.fake_micro_intron_retention,
    }:
        if event.isoform_region != SupplementaryMatchConstants.undefined_region:
            add = ":" + regions_to_str(isoform_introns[event.isoform_region[0] : event.isoform_region[1] + 1])
    elif es in {
        MatchEventSubtype.major_exon_elongation_left,
        MatchEventSubtype.major_exon_elongation_right,
        MatchEventSubtype.exon_elongation_left,
        MatchEventSubtype.exon_elongation_right,
        MatchEventSubtype.alternative_tss_left,
        MatchEventSubtype.alternative_tss_right,
        MatchEventSubtype.alternative_polya_site_left,
        MatchEventSubtype.alternative_polya_site_right,
        MatchEventSubtype.correct_polya_site_right,
        MatchEventSubtype.correct_polya_site_left,
        MatchEventSubtype.internal_polya_left,
        MatchEventSubtype.internal_polya_right,
        MatchEventSubtype.terminal_site_match_left,
        MatchEventSubtype.terminal_site_match_right,
        MatchEventSubtype.terminal_site_match_left_precise,
        MatchEventSubtype.terminal_site_match_right_precise,
    }:
        add = f":{event.event_info}"
    else:
        if event.read_region != SupplementaryMatchConstants.undefined_region and event.read_region[0] >= 0 and event.read_region[1] >= 0:
            add = ":" + regions_to_str(read_introns[event.read_region[0] : event.read_region[1] + 1])
    return match_subtype_to_str(event, strand) + add


def is_matching_assignment(isoform_assignment: ReadAssignment) -> bool:
    if isoform_assignment.assignment_type == ReadAssignmentType.unique:
        return True
    if isoform_assignment.assignment_type.is_unique():
        allowed = {
            MatchEventSubtype.none,
            MatchEventSubtype.fsm,
            MatchEventSubtype.exon_misalignment,
            MatchEventSubtype.intron_shift,
            MatchEventSubtype.terminal_site_match_left,
            MatchEventSubtype.terminal_site_match_right,
            MatchEventSubtype.terminal_site_match_left_precise,
            MatchEventSubtype.terminal_site_match_right_precise,
            MatchEventSubtype.correct_polya_site_right,
            MatchEventSubtype.correct_polya_site_left,
            MatchEventSubtype.exon_elongation_left,
            MatchEventSubtype.exon_elongation_right,
        }
        return all(m.event_type in allowed for m in isoform_assignment.isoform_matches[0].match_subclassifications)
    return False


class PolyAFixer:
    def __init__(self, params) -> None:
        self.params = params

    def correct_read_info(self, read_exons: Sequence[Tuple[int, int]], polya_info: PolyAInfo) -> Tuple[int, int]:
        if len(read_exons) == 1:
            return 0, 0
        pa = self.count_polya_exons(read_exons, polya_info.internal_polya_pos)
        pt = self.count_polyt_exons(read_exons, polya_info.internal_polyt_pos)
        if pt + pa == len(read_exons):
            pt -= 1
            pa -= 1
        return pa, pt

    def count_polya_exons(self, read_exons: Sequence[Tuple[int, int]], internal_polya_pos: int) -> int:
        if internal_polya_pos == -1:
            return 0
        cnt = 0
        for i in range(len(read_exons)):
            exon = read_exons[-i - 1]
            if exon[1] <= internal_polya_pos:
                break
            len_to_polya = internal_polya_pos - exon[0]
            if len_to_polya <= 0 or (
                len_to_polya <= self.params.max_fake_terminal_exon_len and exon[1] - internal_polya_pos > 2 * len_to_polya
            ):
                cnt += 1
        return cnt

    def count_polyt_exons(self, read_exons: Sequence[Tuple[int, int]], internal_polyt_pos: int) -> int:
        if internal_polyt_pos == -1:
            return 0
        cnt = 0
        for i in range(len(read_exons)):
            exon = read_exons[i]
            if exon[0] >= internal_polyt_pos:
                break
            len_to_polyt = exon[1] - internal_polyt_pos
            if len_to_polyt <= 0 or (
                len_to_polyt <= self.params.max_fake_terminal_exon_len and 2 * len_to_polyt < internal_polyt_pos - exon[0]
            ):
                cnt += 1
        return cnt


class PolyAVerifier:
    """Validates read 3' polyA/5' polyT consistency against transcript terminal exons."""

    def __init__(self, gene_info, params) -> None:
        self.params = params
        self.gene_info = gene_info
        self.polya_fixer = PolyAFixer(params)

    def verify_read_ends(self, combined_read_profile, isoform_id: Optional[str], matching_events: List[MatchEvent]):
        if isoform_id is None or combined_read_profile.polya_info is None:
            return matching_events
        read_exons = combined_read_profile.read_split_exon_profile.read_features
        polya_info = combined_read_profile.polya_info
        isoform_exons = self.gene_info.all_isoforms_exons[isoform_id]
        strand = self.gene_info.isoform_strands[isoform_id]

        if strand == "+":
            matching_events, is_internal = self.check_internal_polya(polya_info.internal_polya_pos, matching_events)
            if not is_internal and (polya_info.external_polya_pos != -1 or polya_info.internal_polya_pos != -1):
                matching_events = self.verify_polya(isoform_exons, read_exons, polya_info, matching_events)
        elif strand == "-":
            matching_events, is_internal = self.check_internal_polyt(polya_info.internal_polyt_pos, matching_events)
            if not is_internal and (polya_info.external_polyt_pos != -1 or polya_info.internal_polyt_pos != -1):
                matching_events = self.verify_polyt(isoform_exons, read_exons, polya_info, matching_events)

        return matching_events if matching_events else [MatchEvent(MatchEventSubtype.none)]

    def verify_polya(self, isoform_exons, read_exons, polya_info, matching_events):
        isoform_end = isoform_exons[-1][1]
        event_to_remove = -1
        fake_cnt = 0
        misaligned = 0
        for i, ev in enumerate(matching_events):
            if ev.event_type in (MatchEventSubtype.major_exon_elongation_right, MatchEventSubtype.exon_elongation_right):
                event_to_remove = i
            elif ev.event_type == MatchEventSubtype.fake_terminal_exon_right:
                fake_cnt += 1
            elif ev.event_type == MatchEventSubtype.terminal_exon_misalignment_right:
                misaligned += 1
        if event_to_remove != -1:
            del matching_events[event_to_remove]
        new_events = self.check_if_close(
            isoform_end, polya_info.external_polya_pos, polya_info.internal_polya_pos, matching_events, MatchEventSubtype.correct_polya_site_right
        )
        if new_events is not None:
            return new_events
        ext_pa, int_pa = self.correct_polya_positions(read_exons, fake_cnt, polya_info)
        if misaligned:
            ext_pa = int_pa = isoform_exons[-1][1]
        else:
            matching_events, ext_pa, int_pa = self.detect_reference_exons_beyond_polya(
                isoform_exons, ext_pa, int_pa, matching_events
            )
        new_events = self.check_if_close(
            isoform_end, ext_pa, int_pa, matching_events, MatchEventSubtype.correct_polya_site_right
        )
        if new_events is not None:
            return new_events
        pos = ext_pa if int_pa == -1 else int_pa
        ev_type = (
            MatchEventSubtype.alternative_polya_site_right
            if abs(pos - isoform_end) > self.params.apa_delta
            else MatchEventSubtype.correct_polya_site_right
        )
        matching_events.append(MatchEvent(ev_type, event_info=pos))
        return matching_events

    def verify_polyt(self, isoform_exons, read_exons, polya_info, matching_events):
        isoform_start = isoform_exons[0][0]
        event_to_remove = -1
        fake_cnt = 0
        misaligned = 0
        for i, ev in enumerate(matching_events):
            if ev.event_type in (MatchEventSubtype.major_exon_elongation_left, MatchEventSubtype.exon_elongation_left):
                event_to_remove = i
            elif ev.event_type == MatchEventSubtype.fake_terminal_exon_left:
                fake_cnt += 1
            elif ev.event_type == MatchEventSubtype.terminal_exon_misalignment_left:
                misaligned += 1
        if event_to_remove != -1:
            del matching_events[event_to_remove]
        new_events = self.check_if_close(
            isoform_start, polya_info.external_polyt_pos, polya_info.internal_polyt_pos, matching_events, MatchEventSubtype.correct_polya_site_left
        )
        if new_events is not None:
            return new_events
        ext_pt, int_pt = self.correct_polyt_positions(read_exons, fake_cnt, polya_info)
        if misaligned:
            ext_pt = int_pt = isoform_exons[0][0]
        else:
            matching_events, ext_pt, int_pt = self.detect_reference_exons_before_polyt(
                isoform_exons, ext_pt, int_pt, matching_events
            )
        new_events = self.check_if_close(
            isoform_start, ext_pt, int_pt, matching_events, MatchEventSubtype.correct_polya_site_left
        )
        if new_events is not None:
            return new_events
        pos = ext_pt if int_pt == -1 else int_pt
        ev_type = (
            MatchEventSubtype.alternative_polya_site_left
            if abs(pos - isoform_start) > self.params.apa_delta
            else MatchEventSubtype.correct_polya_site_left
        )
        matching_events.append(MatchEvent(ev_type, event_info=pos))
        return matching_events

    def check_internal_polya(self, internal_polya_pos: int, matching_events: List[MatchEvent]):
        if internal_polya_pos == -1:
            return matching_events, False
        for ev in matching_events:
            if ev.event_type == MatchEventSubtype.incomplete_intron_retention_right:
                matching_events.append(
                    MatchEvent(MatchEventSubtype.internal_polya_right, isoform_region=ev.isoform_region, event_info=internal_polya_pos)
                )
                return matching_events, True
        return matching_events, False

    def check_internal_polyt(self, internal_polyt_pos: int, matching_events: List[MatchEvent]):
        if internal_polyt_pos == -1:
            return matching_events, False
        for ev in matching_events:
            if ev.event_type == MatchEventSubtype.incomplete_intron_retention_left:
                matching_events.append(
                    MatchEvent(MatchEventSubtype.internal_polya_left, isoform_region=ev.isoform_region, event_info=internal_polyt_pos)
                )
                return matching_events, True
        return matching_events, False

    def correct_polya_positions(self, read_exons, fake_cnt: int, polya_info: PolyAInfo) -> Tuple[int, int]:
        assert fake_cnt < len(read_exons)
        return (
            shift_polya(read_exons, fake_cnt, polya_info.external_polya_pos),
            shift_polya(read_exons, fake_cnt, polya_info.internal_polya_pos),
        )

    def correct_polyt_positions(self, read_exons, fake_cnt: int, polya_info: PolyAInfo) -> Tuple[int, int]:
        assert fake_cnt < len(read_exons)
        return (
            shift_polyt(read_exons, fake_cnt, polya_info.external_polyt_pos),
            shift_polyt(read_exons, fake_cnt, polya_info.internal_polyt_pos),
        )

    def detect_reference_exons_beyond_polya(self, isoform_exons, external_polya_pos, internal_polya_pos, matching_events):
        terminal_exon_count = 0
        polya_pos = internal_polya_pos if internal_polya_pos != -1 else external_polya_pos
        while terminal_exon_count < len(isoform_exons) and isoform_exons[-terminal_exon_count - 1][0] >= polya_pos:
            terminal_exon_count += 1
        if terminal_exon_count == len(isoform_exons) or terminal_exon_count == 0:
            return matching_events, external_polya_pos, internal_polya_pos
        term_len = intervals_total_length(isoform_exons[-terminal_exon_count:])
        d_ext = abs(isoform_exons[-terminal_exon_count - 1][1] - external_polya_pos)
        d_int = abs(isoform_exons[-terminal_exon_count - 1][1] - internal_polya_pos)
        dist = min(d_ext, d_int)
        if (term_len <= self.params.max_fake_terminal_exon_len and dist <= self.params.max_fake_terminal_exon_len) or (
            term_len <= self.params.max_missed_exon_len and abs(term_len - dist) <= self.params.delta
        ):
            for i in range(terminal_exon_count):
                rgn = (len(isoform_exons) - 2 - i, len(isoform_exons) - 2 - i)
                matching_events.append(MatchEvent(MatchEventSubtype.terminal_exon_misalignment_right, rgn))
            end_c = isoform_exons[-1][1]
            return matching_events, end_c, end_c
        return matching_events, external_polya_pos, internal_polya_pos

    def detect_reference_exons_before_polyt(self, isoform_exons, external_polyt_pos, internal_polyt_pos, matching_events):
        terminal_exon_count = 0
        polyt_pos = internal_polyt_pos if internal_polyt_pos != -1 else external_polyt_pos
        while terminal_exon_count < len(isoform_exons) and isoform_exons[terminal_exon_count][1] <= polyt_pos:
            terminal_exon_count += 1
        if terminal_exon_count == 0 or terminal_exon_count == len(isoform_exons):
            return matching_events, external_polyt_pos, internal_polyt_pos
        term_len = intervals_total_length(isoform_exons[:terminal_exon_count])
        d_ext = abs(isoform_exons[terminal_exon_count][0] - external_polyt_pos)
        d_int = abs(isoform_exons[terminal_exon_count][0] - internal_polyt_pos)
        dist = min(d_ext, d_int)
        if (term_len <= self.params.max_fake_terminal_exon_len and dist <= self.params.max_fake_terminal_exon_len) or (
            term_len <= self.params.max_missed_exon_len and abs(term_len - dist) <= self.params.delta
        ):
            for i in range(terminal_exon_count):
                matching_events.append(MatchEvent(MatchEventSubtype.terminal_exon_misalignment_left, (i, i)))
            st_c = isoform_exons[0][0]
            return matching_events, st_c, st_c
        return matching_events, external_polyt_pos, internal_polyt_pos

    def check_if_close(self, isoform_bound: int, external_pos: int, internal_pos: int, matching_events, event_type):
        d_ext = abs(isoform_bound - external_pos) if external_pos != -1 else math.inf
        d_int = abs(isoform_bound - internal_pos) if internal_pos != -1 else math.inf
        if d_int <= self.params.apa_delta and d_int <= d_ext:
            matching_events.append(MatchEvent(event_type, event_info=internal_pos))
            return matching_events
        if d_ext <= self.params.apa_delta and d_int > d_ext:
            matching_events.append(MatchEvent(event_type, event_info=external_pos))
            return matching_events
        return None
