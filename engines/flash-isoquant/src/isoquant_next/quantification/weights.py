"""Counting strategy weights, assignment feature extractors, and ungrouped TPM scaling."""

from __future__ import annotations

import logging
from typing import Iterable, Set

from isoquant_next.assignment.types import MatchEventSubtype, ReadAssignmentType
from isoquant_next.config.enums import CountingStrategy, NormalizationMethod

logger = logging.getLogger("IsoQuant")

# Full-length unique match events for fsm_only: spliced FSM or mono-exon read on mono-exon isoform
_FSM_ONLY_EVENTS = frozenset({MatchEventSubtype.fsm, MatchEventSubtype.mono_exon_match})


class CountingStrategyFlags:
    """Boolean feature/weight switches decoded from CountingStrategy."""

    __slots__ = ("use_ambiguous", "use_inconsistent_minor", "use_inconsistent", "fsm_only")

    def __init__(self, counting_strategy: CountingStrategy) -> None:
        self.use_ambiguous: bool = counting_strategy.ambiguous()
        self.use_inconsistent_minor: bool = counting_strategy.inconsistent_minor()
        self.use_inconsistent: bool = counting_strategy.inconsistent()
        self.fsm_only: bool = counting_strategy == CountingStrategy.fsm_only


class GeneAssignmentExtractor:
    """Extract gene-level feature IDs and assignment types from a ReadAssignment."""

    @staticmethod
    def get_features(read_assignment) -> Set[str]:
        return {m.assigned_gene for m in read_assignment.isoform_matches if m.assigned_gene}

    @staticmethod
    def get_assignment_type(read_assignment):
        return read_assignment.gene_assignment_type

    @staticmethod
    def confirms_feature(read_assignment) -> bool:
        return read_assignment.gene_assignment_type.is_unique()


class TranscriptAssignmentExtractor:
    """Extract transcript-level feature IDs and assignment types from a ReadAssignment."""

    @staticmethod
    def get_features(read_assignment) -> Set[str]:
        return {m.assigned_transcript for m in read_assignment.isoform_matches if m.assigned_transcript}

    @staticmethod
    def get_assignment_type(read_assignment):
        return read_assignment.assignment_type

    @staticmethod
    def confirms_feature(read_assignment) -> bool:
        transcript_id = read_assignment.isoform_matches[0].assigned_transcript
        return (
            read_assignment.assignment_type.is_unique()
            and (
                len(read_assignment.gene_info.all_isoforms_introns[transcript_id]) == 0
                or len(read_assignment.corrected_exons) > 1
            )
        )


class ReadWeightCounter:
    """Compute fractional read weights according to the selected CountingStrategy."""

    def __init__(self, strategy: str | CountingStrategy) -> None:
        self.strategy: CountingStrategy = (
            strategy if isinstance(strategy, CountingStrategy) else CountingStrategy[strategy]
        )
        self.strategy_flags = CountingStrategyFlags(self.strategy)

    def process_inconsistent(self, assignment_type: ReadAssignmentType, feature_count: int) -> float:
        if assignment_type == ReadAssignmentType.inconsistent_ambiguous or feature_count > 1:
            if self.strategy_flags.use_ambiguous and self.strategy_flags.use_inconsistent:
                return 1.0 / feature_count
            return 0.0
        if self.strategy_flags.use_inconsistent:
            return 1.0
        if (
            self.strategy_flags.use_inconsistent_minor
            and assignment_type == ReadAssignmentType.inconsistent_non_intronic
        ):
            return 1.0
        return 0.0

    def process_ambiguous(self, feature_count: int) -> float:
        if feature_count == 0:
            return 0.0
        if feature_count == 1:
            return 1.0
        if self.strategy_flags.use_ambiguous:
            return 1.0 / float(feature_count)
        return 0.0

    def process_unique(self, read_assignment) -> float:
        if not self.strategy_flags.fsm_only:
            return 1.0
        is_full_match = any(
            e.event_type in _FSM_ONLY_EVENTS
            for m in read_assignment.isoform_matches
            for e in m.match_subclassifications
        )
        return 1.0 if is_full_match else 0.0


def convert_ungrouped_to_tpm(
    counts_file_name: str,
    output_tpm_file_name: str,
    normalization: NormalizationMethod = NormalizationMethod.simple,
    reads_for_tpm: int = -1,
) -> None:
    """Convert ungrouped count TSV to TPM TSV without loading the table into DataFrames."""
    total_counts = 0.0
    with open(counts_file_name, "r", encoding="utf-8") as inf:
        for line in inf:
            if line.startswith("_") or line.startswith("feature_id\t"):
                continue
            fs = line.rstrip("\r\n").split("\t")
            if len(fs) >= 2:
                total_counts += float(fs[1])

    if total_counts == 0.0:
        logger.error("No counts were detected in %s, cannot convert to TPM", counts_file_name)
        return

    unassigned_tpm = 0.0
    if normalization == NormalizationMethod.usable_reads and reads_for_tpm > 0:
        total_reads = float(reads_for_tpm)
        unassigned_tpm = 1000000.0 * (1.0 - total_counts / total_reads)
    else:
        total_reads = total_counts
    scale_factor = 1000000.0 / total_reads

    with open(output_tpm_file_name, "w", encoding="utf-8") as outf, open(
        counts_file_name, "r", encoding="utf-8"
    ) as inf:
        for line in inf:
            if line.startswith("_"):
                break
            if line.startswith("feature_id\t"):
                outf.write(line.replace("count", "TPM"))
                continue
            fs = line.rstrip("\r\n").split("\t")
            if len(fs) >= 2:
                feature_id, count = fs[0], float(fs[1])
                outf.write(f"{feature_id}\t{scale_factor * count:.6f}\n")
        if unassigned_tpm > 0.0:
            outf.write(f"__unassigned\t{unassigned_tpm:.6f}\n")


__all__ = [
    "CountingStrategyFlags",
    "GeneAssignmentExtractor",
    "ReadWeightCounter",
    "TranscriptAssignmentExtractor",
    "convert_ungrouped_to_tpm",
]
