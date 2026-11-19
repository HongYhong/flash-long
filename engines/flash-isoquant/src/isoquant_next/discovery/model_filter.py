"""Stage 4a of transcript discovery: pre-filtering, structural filtering, and splice-site correction."""

from __future__ import annotations

import logging
from typing import Dict, List

from isoquant_next.alignment.alignment_info import PolyAInfo
from isoquant_next.annotation.gene_info import GeneInfo, TranscriptModelType
from isoquant_next.assignment.assigner import LongReadAssigner
from isoquant_next.assignment.profiles import CombinedProfileConstructor
from isoquant_next.assignment.types import ReadAssignmentType, is_matching_assignment
from isoquant_next.discovery.intron_graph import TerminalVertex
from isoquant_next.discovery.splice_corrector import (
    SUPPORTED_STRANDS,
    correct_splice_site_errors,
    count_deletions_for_splice_site_locations,
    generate_updated_exon_list,
)

logger = logging.getLogger("IsoQuant")


class ModelFilter:
    """Structural filtering & splice-site correction of discovered transcript models."""

    def __init__(self, ctx, end_processor) -> None:
        self.ctx = ctx
        self.store = ctx.store
        self.args = ctx.args
        self.chr_record = ctx.chr_record
        self.string_pools = ctx.string_pools
        self.create_nics = ctx.create_alt_ts_nics
        self.novel_apa = ctx.novel_apa
        self.use_tss_model = ctx.use_tss_model
        self.end_processor = end_processor

    def pre_filter_transcripts(self) -> None:
        counts = sorted(
            (
                self.store.internal_counter[m.transcript_id]
                for m in self.store.transcript_model_storage
                if len(m.exon_blocks) <= 2
            ),
            reverse=True,
        )
        coverage_cutoff = counts[50] if len(counts) > 50 else self.args.min_novel_count

        filtered = []
        for model in self.store.transcript_model_storage:
            if len(model.exon_blocks) > 2:
                filtered.append(model)
                continue
            if model.transcript_type != TranscriptModelType.known and (
                self.store.internal_counter[model.transcript_id] < coverage_cutoff
                or self.store.mapping_quality(model.transcript_id) < self.args.simple_models_mapq_cutoff
            ):
                self.store.delete_from_storage(model.transcript_id)
                continue
            filtered.append(model)
        self.store.transcript_model_storage = filtered

    def filter_transcripts(self) -> None:
        pre_filtered = []
        to_substitute = self.detect_similar_isoforms(self.store.transcript_model_storage)

        for model in self.store.transcript_model_storage:
            if model.transcript_type == TranscriptModelType.known:
                pre_filtered.append(model)
                continue

            comp_cov = self.ctx.intron_graph.get_max_component_coverage(model.intron_path)
            if (
                comp_cov == 0
                or not model.intron_path
                or (len(model.intron_path) == 1 and self.ctx.intron_graph.is_monointron(model.intron_path[0]))
            ):
                comp_cov = self.ctx.intron_graph.get_overlapping_component_max_coverage(
                    (model.get_start(), model.get_end())
                )
                cutoff = max(self.args.min_novel_count, self.args.min_mono_count_rel * comp_cov)
            else:
                cutoff = max(self.args.min_novel_count, self.args.min_novel_count_rel * comp_cov)

            if (
                model.transcript_id in to_substitute
                or self.store.internal_counter[model.transcript_id] < cutoff
                or (
                    len(model.exon_blocks) <= 2
                    and self.store.mapping_quality(model.transcript_id) < self.args.simple_models_mapq_cutoff
                )
            ):
                self.store.delete_from_storage(model.transcript_id)
                continue

            self.end_processor.correct_novel_transcript_ends(
                model, self.store.transcript_read_ids[model.transcript_id]
            )
            self.end_processor._mark_terminal_confirmation(model)
            pre_filtered.append(model)

        to_substitute = self.detect_similar_isoforms(pre_filtered, remove_unconfirmed=True)
        filtered = []
        for model in pre_filtered:
            if model.transcript_type != TranscriptModelType.known and model.transcript_id in to_substitute:
                self.store.delete_from_storage(model.transcript_id)
                continue
            filtered.append(model)

        self.store.transcript_model_storage = filtered
        if self.create_nics or self.novel_apa:
            self.end_processor._add_known_alternative_end_models()
            self.end_processor._drop_duplicate_alt_end_models()

    def correct_transcript_splice_sites(self) -> None:
        if getattr(self.args, "no_splice_site_correction", False):
            return
        for model in self.store.transcript_model_storage:
            if model.transcript_type == TranscriptModelType.known:
                continue
            corrected = self._correct_one_model_splice_sites(
                model, self.store.transcript_read_ids[model.transcript_id]
            )
            if corrected and corrected != model.exon_blocks and self._valid_exon_chain(corrected):
                model.exon_blocks = corrected

    def _correct_one_model_splice_sites(self, model, reads):
        exons = model.exon_blocks
        if len(exons) < 2 or not reads or model.strand not in SUPPORTED_STRANDS:
            return None
        splice_cases = {}
        for read in reads:
            if not read.exons:
                continue
            count_deletions_for_splice_site_locations(
                read.exons[0][0], read.exons[-1][1], read.boundary_deletions or [], exons, splice_cases
            )
        if not splice_cases:
            return None
        locs = correct_splice_site_errors(splice_cases, self.chr_record, model.strand)
        if not locs:
            return None
        return generate_updated_exon_list(splice_cases, locs, exons)

    @staticmethod
    def _valid_exon_chain(exons) -> bool:
        for i, (s, e) in enumerate(exons):
            if s > e or (i > 0 and s <= exons[i - 1][1]):
                return False
        return True

    def detect_similar_isoforms(self, model_storage, remove_unconfirmed: bool = False) -> Dict[str, str]:
        to_substitute: Dict[str, str] = {}
        for model in model_storage:
            if len(model.exon_blocks) <= 2 or model.transcript_id in to_substitute:
                continue
            gi = GeneInfo.from_models([model], self.args.delta)
            assigner = LongReadAssigner(gi, self.args, self.string_pools)
            pc = CombinedProfileConstructor(gi, self.args)

            for m in model_storage:
                if (
                    m.transcript_type == TranscriptModelType.known
                    or m.transcript_id == model.transcript_id
                    or m.transcript_id in to_substitute
                    or len(m.exon_blocks) == 1
                    or not m.intron_path
                    or len(m.exon_blocks) > len(model.exon_blocks)
                ):
                    continue

                if m.intron_path[0][0] == TerminalVertex.polyt:
                    pa = PolyAInfo(-1, m.intron_path[0][1], -1, -1)
                elif m.intron_path[-1][0] == TerminalVertex.polya:
                    pa = PolyAInfo(m.intron_path[-1][1], -1, -1, -1)
                else:
                    pa = PolyAInfo(-1, -1, -1, -1)
                cp = pc.construct_profiles(m.exon_blocks, pa, [])
                ra = assigner.assign_to_isoform(m.transcript_id, cp)
                if is_matching_assignment(ra):
                    to_substitute[m.transcript_id] = model.transcript_id
                elif (
                    remove_unconfirmed
                    and ra.assignment_type
                    in (ReadAssignmentType.unique_minor_difference, ReadAssignmentType.inconsistent_non_intronic)
                    and not (m.polya_confirmed or (self.use_tss_model and m.tss_confirmed))
                ):
                    to_substitute[m.transcript_id] = model.transcript_id

        return to_substitute
