"""Multimap assignment resolution supporting CompactMultimapTuple and BasicReadAssignment."""

from __future__ import annotations

import logging
import math
from typing import Iterable, List, Optional, Sequence, Tuple, TypeVar

from isoquant_next.assignment.types import CompactMultimapTuple, ReadAssignmentType
from isoquant_next.config.enums import MultimapResolvingStrategy
from isoquant_next.core.utils import intersection_len

logger = logging.getLogger("IsoQuant")

TAssignment = TypeVar("TAssignment")


def _same_assignment(a, b) -> bool:
    """Compare two multimap records without hydrating string pools when both are compact tuples."""
    if isinstance(a, CompactMultimapTuple) and isinstance(b, CompactMultimapTuple):
        return (
            a.read_id == b.read_id
            and a.chr_id_int == b.chr_id_int
            and a.start == b.start
            and a.end == b.end
            and a.isoform_ids == b.isoform_ids
        )
    return a == b


def _extract_feature_ids(item) -> Tuple[Sequence, Sequence]:
    """Return (gene_keys, isoform_keys) preferring integer IDs over Python strings."""
    gene_ids = getattr(item, "gene_ids", None)
    isoform_ids = getattr(item, "isoform_ids", None)
    if gene_ids or isoform_ids or isinstance(item, CompactMultimapTuple):
        return gene_ids or (), isoform_ids or ()
    return getattr(item, "genes", ()), getattr(item, "isoforms", ())


class MultimapResolver:
    """Resolves secondary/multimap read assignments without cross-process object pickling."""

    duplicate_counter = 0

    def __init__(self, strategy: MultimapResolvingStrategy) -> None:
        self.strategy = strategy

    def resolve(self, assignment_list: Optional[List[TAssignment]]) -> Optional[List[TAssignment]]:
        if assignment_list is None or len(assignment_list) <= 1:
            return assignment_list

        if self.strategy == MultimapResolvingStrategy.ignore_multimapper:
            for a in assignment_list:
                a.assignment_type = ReadAssignmentType.suspended
            return assignment_list
        if self.strategy == MultimapResolvingStrategy.merge:
            return self.merge_assignments(assignment_list)
        if self.strategy == MultimapResolvingStrategy.take_best:
            return self.select_best_assignment(assignment_list)
        raise ValueError("Unsupported multimap strategy")

    def select_best_assignment(self, assignment_list: List[TAssignment]) -> List[TAssignment]:
        logger.debug("Resolving read %s", assignment_list[0].read_id)
        primary_unique: List[int] = []
        consistent_assignments: List[int] = []
        inconsistent_assignments: List[int] = []
        primary_inconsistent: List[int] = []
        noninformative: List[int] = []

        for i, a in enumerate(assignment_list):
            if a.assignment_type.is_inconsistent():
                if not a.multimapper:
                    primary_inconsistent.append(i)
                inconsistent_assignments.append(i)
            elif a.assignment_type.is_consistent():
                consistent_assignments.append(i)
                if not a.multimapper and a.assignment_type != ReadAssignmentType.ambiguous:
                    primary_unique.append(i)
            else:
                noninformative.append(i)

        if primary_unique:
            return self.filter_assignments(assignment_list, primary_unique)
        if consistent_assignments:
            return self.filter_assignments(assignment_list, consistent_assignments)
        if primary_inconsistent:
            return self.select_best_inconsistent(assignment_list, primary_inconsistent)
        if inconsistent_assignments:
            return self.select_best_inconsistent(assignment_list, inconsistent_assignments)
        if noninformative:
            return self.select_noninformative(assignment_list, noninformative)

        logger.warning("Unexpected assignments in multimap resolution %s", assignment_list[0].read_id)
        return assignment_list[0:1]

    @staticmethod
    def select_best_inconsistent(
        assignment_list: List[TAssignment], inconsistent_assignments: Sequence[int]
    ) -> List[TAssignment]:
        if len(inconsistent_assignments) <= 1:
            return MultimapResolver.filter_assignments(assignment_list, inconsistent_assignments)

        scores = [(assignment_list[i].penalty_score, i) for i in inconsistent_assignments]
        best_score = min(scores, key=lambda x: x[0])[0]
        best_indices = [idx for score, idx in scores if score == best_score]
        return MultimapResolver.filter_assignments(assignment_list, best_indices)

    @staticmethod
    def merge_assignments(assignment_list: List[TAssignment]) -> List[TAssignment]:
        informative = [
            i
            for i in range(len(assignment_list))
            if assignment_list[i].assignment_type != ReadAssignmentType.noninformative
        ]
        return MultimapResolver.filter_assignments(assignment_list, informative)

    @staticmethod
    def filter_assignments(
        assignment_list: List[TAssignment], assignments_to_keep: Iterable[int]
    ) -> List[TAssignment]:
        kept_list = MultimapResolver.find_duplicates(assignment_list, list(assignments_to_keep))

        all_genes = set()
        all_isoforms = set()
        for i in kept_list:
            g_keys, t_keys = _extract_feature_ids(assignment_list[i])
            all_genes.update(g_keys)
            all_isoforms.update(t_keys)

        change_tx_type = len(all_isoforms) > 1
        change_gene_type = len(all_genes) > 1
        kept_set = set(kept_list)

        for i, assignment in enumerate(assignment_list):
            if i in kept_set:
                amb_type = (
                    ReadAssignmentType.inconsistent_ambiguous
                    if assignment.assignment_type.is_inconsistent()
                    else ReadAssignmentType.ambiguous
                )
                if change_tx_type:
                    assignment.assignment_type = amb_type
                    assignment.multimapper = True
                if change_gene_type:
                    assignment.gene_assignment_type = amb_type
                    assignment.multimapper = True
            else:
                assignment.assignment_type = ReadAssignmentType.suspended
                assignment.gene_assignment_type = ReadAssignmentType.suspended

        return assignment_list

    @staticmethod
    def find_duplicates(
        assignment_list: List[TAssignment], assignment_indices: Iterable[int]
    ) -> List[int]:
        indices = list(assignment_indices)
        if len(indices) <= 1:
            return indices

        selected: List[int] = []
        discarded: set[int] = set()
        for i in range(len(indices)):
            idx1 = indices[i]
            if idx1 in discarded:
                continue
            selected.append(idx1)
            for j in range(i + 1, len(indices)):
                idx2 = indices[j]
                if idx2 in discarded:
                    continue
                if _same_assignment(assignment_list[idx1], assignment_list[idx2]):
                    discarded.add(idx2)
                    if assignment_list[idx1].genomic_region == assignment_list[idx2].genomic_region:
                        MultimapResolver.duplicate_counter += 1
                        ex = assignment_list[idx1]
                        chr_label = ex.chr_id if hasattr(ex, "chr_id") else str(getattr(ex, "chr_id_int", "."))
                        if MultimapResolver.duplicate_counter == 5:
                            logger.warning(
                                "More possible duplicates were detected but will not be printed. "
                                "Use --debug to see all of them in the isoquant.log. "
                                "All duplicated entries will be ignored."
                            )
                        msg = (
                            f"Read {ex.read_id} seems to have duplicated alignment entries at "
                            f"{chr_label}: {ex.start} and will be ignored."
                        )
                        if MultimapResolver.duplicate_counter >= 5:
                            logger.debug(msg)
                        else:
                            logger.info(msg)
        return selected

    @staticmethod
    def select_noninformative(
        assignment_list: List[TAssignment], assignment_indices: Sequence[int]
    ) -> List[TAssignment]:
        overlap_index_list = []
        max_overlap_len = 0
        for i in assignment_indices:
            a = assignment_list[i]
            ov = intersection_len(a.genomic_region, (a.start, a.end))
            if ov > max_overlap_len:
                max_overlap_len = ov
            overlap_index_list.append((ov, a.genomic_region[0], i))

        min_region_start = math.inf
        best_assignment = -1
        for ov, reg_start, idx in overlap_index_list:
            if ov == max_overlap_len and reg_start < min_region_start:
                min_region_start = reg_start
                best_assignment = idx

        assert best_assignment != -1
        return MultimapResolver.filter_assignments(assignment_list, [best_assignment])


def resolve_compact_multimappers(
    records: List[CompactMultimapTuple], strategy: MultimapResolvingStrategy
) -> List[CompactMultimapTuple]:
    """Zero-pickle resolution over fixed-width CompactMultimapTuple records."""
    resolved = MultimapResolver(strategy).resolve(records)
    return resolved if resolved is not None else []


def resolve_compact_bytes(payloads: Sequence[bytes], strategy: MultimapResolvingStrategy) -> List[bytes]:
    """Resolve binary-packed CompactMultimapTuple payloads without hydrating ReadAssignment."""
    records = [CompactMultimapTuple.from_bytes(b) for b in payloads]
    resolve_compact_multimappers(records, strategy)
    return [r.to_bytes() for r in records]
