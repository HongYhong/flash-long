"""Stages 4b & 4c of transcript discovery: gated read assignment, ownership resolution, and counting."""

from __future__ import annotations

import copy
import logging
from collections import defaultdict
from typing import Dict, List

from isoquant_next.annotation.gene_info import GeneInfo, TranscriptModelType
from isoquant_next.assignment.assigner import LongReadAssigner
from isoquant_next.assignment.profiles import CombinedProfileConstructor
from isoquant_next.assignment.types import (
    IsoformMatch,
    MatchClassification,
    ReadAssignment,
    ReadAssignmentType,
)
from isoquant_next.discovery.state import reference_unique_known_isoform

logger = logging.getLogger("IsoQuant")


class ConstructionReadAssigner:
    """Reassigns reads to surviving models and resolves FL-vs-nonFL ownership conflicts."""

    def __init__(self, ctx) -> None:
        self.ctx = ctx
        self.store = ctx.store
        self.args = ctx.args
        self.string_pools = ctx.string_pools

    def resolve_read_ownership_conflicts(self) -> None:
        if not self.store.nonfl_known_reads:
            return
        owned: Dict[str, int] = defaultdict(int)
        for mid in self.store.construction_assignment.values():
            owned[mid] += 1

        read_to_fl: Dict[str, str] = {}
        for model in self.store.transcript_model_storage:
            if model.transcript_id in self.store.nonfl_known_reads:
                continue
            if model.transcript_type == TranscriptModelType.known and owned.get(model.transcript_id, 0) > 0:
                continue
            for a in self.store.transcript_read_ids[model.transcript_id]:
                read_to_fl.setdefault(a.read_id, model.transcript_id)

        surviving = {m.transcript_id for m in self.store.transcript_model_storage}
        emptied: List[str] = []
        for known_id, reads in self.store.nonfl_known_reads.items():
            if known_id not in surviving:
                continue
            for a in reads:
                if self.store.construction_assignment.get(a.read_id) != known_id:
                    continue
                fl_id = read_to_fl.get(a.read_id)
                if fl_id is None or fl_id == known_id:
                    continue
                self.store.construction_assignment[a.read_id] = fl_id
                if a in self.store.transcript_read_ids[known_id]:
                    self.store.transcript_read_ids[known_id].remove(a)
                    self.store.internal_counter[known_id] -= 1
                    self.store.read_assignment_counts[a.read_id] -= 1
            if self.store.internal_counter[known_id] <= 0:
                emptied.append(known_id)

        if emptied:
            drop = set(emptied)
            for kid in drop:
                self.store.delete_from_storage(kid)
            self.store.transcript_model_storage = [
                m for m in self.store.transcript_model_storage if m.transcript_id not in drop
            ]

    def assign_reads_to_models(self, read_assignments) -> None:
        if not self.store.transcript_model_storage:
            for a in read_assignments:
                self.store.read_assignment_counts[a.read_id] = 0
            return

        model_gi = GeneInfo.from_models(self.store.transcript_model_storage, self.args.delta)
        assigner = LongReadAssigner(model_gi, self.args, self.string_pools, quick_mode=True)
        pc = CombinedProfileConstructor(model_gi, self.args)
        model_ids = set(model_gi.all_isoforms_exons.keys())

        for a in read_assignments:
            read_id = a.read_id
            if self.store.read_assignment_counts[read_id] > 0:
                continue

            known_iso = reference_unique_known_isoform(a)
            if known_iso is not None:
                if known_iso in model_ids:
                    self.store.transcript_read_ids[known_iso].append(a)
                    self.store.internal_counter[known_iso] += 1
                    self.store.read_assignment_counts[read_id] += 1
                else:
                    self.store.read_assignment_counts[read_id] = 0
                continue

            cp = pc.construct_profiles(a.corrected_exons, a.polya_info, [])
            ma = assigner.assign_to_isoform(read_id, cp)
            ma.read_group = a.read_group
            if ma.assignment_type.is_consistent():
                matched = [m.assigned_transcript for m in ma.isoform_matches]
                if len(matched) == 1:
                    self.store.internal_counter[matched[0]] += 1
                for m in ma.isoform_matches:
                    self.store.read_assignment_counts[read_id] += 1
                    self.store.transcript_read_ids[m.assigned_transcript].append(a)
            else:
                self.store.read_assignment_counts[read_id] = 0


class ModelReadCounter:
    """Produces final per-read ReadAssignments against the discovered model set & feeds counters."""

    def __init__(self, ctx, transcript_counter, gene_counter) -> None:
        self.ctx = ctx
        self.store = ctx.store
        self.gene_info = ctx.gene_info
        self.args = ctx.args
        self.string_pools = ctx.string_pools
        self.transcript_counter = transcript_counter
        self.gene_counter = gene_counter
        self.model_read_assignments: List[ReadAssignment] = []
        self._graph_tolerance_args = None

    def _single_model_assigner(self, model, cache):
        entry = cache.get(model.transcript_id)
        if entry is None:
            gi = GeneInfo.from_models([model], self.args.delta)
            entry = (
                LongReadAssigner(gi, self._graph_tolerance_args, self.string_pools, quick_mode=False),
                CombinedProfileConstructor(gi, self._graph_tolerance_args),
            )
            cache[model.transcript_id] = entry
        return entry

    def _model_read_assignment_single(self, source, transcript_id: str, model_by_id, cache) -> ReadAssignment:
        model = model_by_id[transcript_id]
        if model.transcript_type == TranscriptModelType.known and source.assignment_type in (
            ReadAssignmentType.unique,
            ReadAssignmentType.unique_minor_difference,
        ):
            reuse = next((m for m in source.isoform_matches if m.assigned_transcript == transcript_id), None)
            if reuse is not None:
                return ReadAssignment(source.read_id, source.assignment_type, self.string_pools, match=[reuse])

        assigner, pc = self._single_model_assigner(model, cache)
        profile = pc.construct_profiles(source.corrected_exons, source.polya_info, [])
        ra = assigner.assign_to_isoform(source.read_id, profile)
        if (
            ra.assignment_type.is_consistent()
            and ra.isoform_matches
            and ra.isoform_matches[0].assigned_transcript == transcript_id
        ):
            return ra
        forced = IsoformMatch(
            MatchClassification.full_splice_match,
            self.string_pools,
            assigned_gene=model.gene_id,
            assigned_transcript=transcript_id,
        )
        return ReadAssignment(source.read_id, ReadAssignmentType.unique_minor_difference, self.string_pools, match=[forced])

    def _build_discarded_assignment(self, read_id: str, gene: Optional[str], surviving_genes: set[str]) -> ReadAssignment:
        if gene is not None and gene in surviving_genes:
            match = IsoformMatch(MatchClassification.genic, self.string_pools, assigned_gene=gene)
            ra = ReadAssignment(read_id, ReadAssignmentType.discarded, self.string_pools, match=[match])
            ra.gene_assignment_type = ReadAssignmentType.unique
            return ra
        return ReadAssignment(read_id, ReadAssignmentType.discarded, self.string_pools)

    def build_model_read_assignments(self, read_assignments) -> None:
        self.model_read_assignments = []
        self._graph_tolerance_args = copy.copy(self.args)
        gcd = getattr(self.args, "graph_clustering_distance", None)
        self._graph_tolerance_args.delta = max(self.args.delta, gcd) if gcd else self.args.delta
        tx2gene = {t.transcript_id: t.gene_id for t in self.store.transcript_model_storage}
        model_gi = GeneInfo.from_models(self.store.transcript_model_storage, self.args.delta)
        model_by_id = {t.transcript_id: t for t in self.store.transcript_model_storage}
        model_ids = set(model_by_id.keys())
        surviving_genes = set(tx2gene.values())

        single_cache = {}
        full_assigner = None
        full_pc = None

        for source in read_assignments:
            read_id = source.read_id
            constr_tid = self.store.construction_assignment.get(read_id)
            if constr_tid is not None and constr_tid in model_ids:
                ra = self._model_read_assignment_single(source, constr_tid, model_by_id, single_cache)
            else:
                known_iso = reference_unique_known_isoform(source)
                if known_iso is not None and known_iso in model_ids:
                    ra = self._model_read_assignment_single(source, known_iso, model_by_id, single_cache)
                elif known_iso is not None:
                    gene = self.gene_info.gene_id_map.get(known_iso)
                    ra = self._build_discarded_assignment(read_id, gene, surviving_genes)
                elif model_ids:
                    if full_assigner is None:
                        full_assigner = LongReadAssigner(model_gi, self.args, self.string_pools, quick_mode=False)
                        full_pc = CombinedProfileConstructor(model_gi, self.args)
                    profile = full_pc.construct_profiles(source.corrected_exons, source.polya_info, [])
                    ra = full_assigner.assign_to_isoform(read_id, profile)
                else:
                    ra = ReadAssignment(read_id, ReadAssignmentType.noninformative, self.string_pools)

            ra.gene_info = model_gi
            self._copy_read_aux(ra, source)
            self.model_read_assignments.append(ra)
            if self.transcript_counter is not None:
                self.transcript_counter.add_read_info(ra)
            if self.gene_counter is not None:
                self.gene_counter.add_read_info(ra)

        if self.transcript_counter is not None:
            self.transcript_counter.add_confirmed_features([m.transcript_id for m in self.store.transcript_model_storage])
        if self.gene_counter is not None:
            self.gene_counter.add_confirmed_features({tx2gene[m.transcript_id] for m in self.store.transcript_model_storage})

    @staticmethod
    def _copy_read_aux(ra: ReadAssignment, source: ReadAssignment) -> None:
        ra.chr_id_int = source.chr_id_int
        ra.strand = source.strand
        ra.exons = source.exons
        ra.corrected_exons = source.corrected_exons
        ra.polyA_found = source.polyA_found
        ra.cage_found = source.cage_found
        ra.polya_info = source.polya_info
        ra.barcode_id = source.barcode_id
        ra.umi_id = source.umi_id
        ra.read_group_ids = source.read_group_ids
        ra.additional_attributes = source.additional_attributes
