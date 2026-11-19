"""5-stage GraphBasedModelConstructor orchestrating Node 6 transcript discovery."""

from __future__ import annotations

import copy
import logging
from collections import defaultdict
from typing import Optional

from isoquant_next.alignment.alignment_info import PolyAInfo
from isoquant_next.annotation.gene_info import TranscriptModelType
from isoquant_next.assignment.assigner import LongReadAssigner
from isoquant_next.assignment.profiles import CombinedProfileConstructor
from isoquant_next.assignment.types import (
    IsoformMatch,
    MatchClassification,
    ReadAssignment,
    ReadAssignmentType,
    match_subtype_to_str_with_additional_info,
)
from isoquant_next.core.utils import junctions_from_blocks
from isoquant_next.discovery.assignment_constructor import AssignmentBasedConstructor
from isoquant_next.discovery.end_processor import TranscriptEndProcessor
from isoquant_next.discovery.fl_constructor import GraphBasedConstructor
from isoquant_next.discovery.gene_joiner import TranscriptToGeneJoiner
from isoquant_next.discovery.model_filter import ModelFilter
from isoquant_next.discovery.read_assigner import ConstructionReadAssigner, ModelReadCounter
from isoquant_next.discovery.state import (
    ChromosomeDiscoveryState,
    ModelConstructionContext,
    ModelStore,
    StrandnessReportingLevel,
)

logger = logging.getLogger("IsoQuant")


class GraphBasedModelConstructor:
    """Orchestrates the 5 transcript-discovery stages with explicit ChromosomeDiscoveryState."""

    def __init__(
        self,
        gene_info,
        chr_record,
        args,
        transcript_counter,
        gene_counter,
        id_distributor,
        grouping_strategy_names=None,
        use_technical_replicas: bool = False,
        string_pools=None,
        polya_predictions=None,
        tss_predictions=None,
        chr_state: Optional[ChromosomeDiscoveryState] = None,
    ) -> None:
        store = ModelStore(gene_info, chr_state=chr_state)
        self.ctx = ModelConstructionContext(
            gene_info,
            chr_record,
            args,
            id_distributor,
            string_pools,
            polya_predictions,
            tss_predictions,
            store,
            grouping_strategy_names=grouping_strategy_names,
            use_technical_replicas=use_technical_replicas,
        )
        ctx = self.ctx
        self.store = store
        self.gene_info = ctx.gene_info
        self.args = ctx.args
        self.string_pools = ctx.string_pools
        self.id_distributor = ctx.id_distributor
        self.transcript2transcript = []

        self.model_comparison_assigner, self.model_comparison_profile_constructor = self._make_model_comparators()
        self.fl_constructor = GraphBasedConstructor(ctx)
        self.assignment_constructor = AssignmentBasedConstructor(ctx)
        self.end_processor = TranscriptEndProcessor(ctx)
        self.model_filter = ModelFilter(ctx, self.end_processor)
        self.read_assigner = ConstructionReadAssigner(ctx)
        self.model_read_counter = ModelReadCounter(ctx, transcript_counter, gene_counter)

    def _make_model_comparators(self):
        if not getattr(self.args, "sqanti_output", False):
            return None, None
        exact_args = copy.copy(self.args)
        exact_args.delta = 0
        return (
            LongReadAssigner(self.gene_info, exact_args, self.string_pools),
            CombinedProfileConstructor(self.gene_info, exact_args),
        )

    @property
    def transcript_model_storage(self):
        return self.store.transcript_model_storage

    @property
    def model_read_assignments(self):
        return self.model_read_counter.model_read_assignments

    def get_transcript_id(self):
        return self.id_distributor.increment()

    def process(self, read_assignment_storage) -> None:
        self.ctx.build_graph(read_assignment_storage)
        self.fl_constructor.construct_fl_isoforms()
        self.assignment_constructor.construct_assignment_based_isoforms(read_assignment_storage)
        self.model_filter.pre_filter_transcripts()
        self.read_assigner.assign_reads_to_models(read_assignment_storage)
        self.model_filter.filter_transcripts()
        self.model_filter.correct_transcript_splice_sites()
        self.read_assigner.assign_reads_to_models(read_assignment_storage)
        self.read_assigner.resolve_read_ownership_conflicts()

        joiner = TranscriptToGeneJoiner(self.store.transcript_model_storage, self.gene_info)
        self.store.transcript_model_storage = joiner.join_transcripts()
        self.model_read_counter.build_model_read_assignments(read_assignment_storage)

        if getattr(self.args, "sqanti_output", False):
            self.compare_models_with_known()

    def compare_models_with_known(self) -> None:
        if not self.gene_info.all_isoforms_exons:
            for model in self.store.transcript_model_storage:
                ra = ReadAssignment(
                    model.transcript_id,
                    ReadAssignmentType.intergenic,
                    self.string_pools,
                    match=IsoformMatch(MatchClassification.intergenic, string_pools=self.string_pools),
                )
                pa = (
                    PolyAInfo(-1, model.exon_blocks[0][0], -1, -1)
                    if model.strand == "-"
                    else PolyAInfo(model.exon_blocks[-1][1], -1, -1, -1)
                )
                ra.polya_info = pa
                ra.cage_found = False
                ra.exons = model.exon_blocks
                ra.strand = model.strand
                ra.chr_id = model.chr_id
                ra.set_additional_info("indel_count", "NA")
                ra.set_additional_info("junctions_with_indels", "NA")
                ra.introns_match = False
                ra.gene_info = self.gene_info
                ra.set_additional_info("FSM_class", "C")
                self.transcript2transcript.append(ra)
            return

        gene_to_model = defaultdict(list)
        for model in self.store.transcript_model_storage:
            gene_to_model[model.gene_id].append(model.transcript_id)

        self.transcript2transcript = []
        for model in self.store.transcript_model_storage:
            if model.transcript_type == TranscriptModelType.known:
                continue
            pa = (
                PolyAInfo(-1, model.exon_blocks[0][0], -1, -1)
                if model.strand == "-"
                else PolyAInfo(model.exon_blocks[-1][1], -1, -1, -1)
            )
            cp = self.model_comparison_profile_constructor.construct_profiles(model.exon_blocks, pa, [])
            ra = self.model_comparison_assigner.assign_to_isoform(model.transcript_id, cp)
            if ra is None:
                continue
            ra.polya_info = pa
            ra.cage_found = False
            ra.exons = model.exon_blocks
            ra.strand = model.strand
            ra.chr_id = model.chr_id
            ra.set_additional_info("indel_count", "NA")
            ra.set_additional_info("junctions_with_indels", "NA")
            ra.introns_match = all(e == 1 for e in cp.read_intron_profile.read_profile)
            ra.gene_info = self.gene_info

            if ra.assignment_type.is_unassigned() or not ra.isoform_matches:
                ra.assignment_type = ReadAssignmentType.intergenic
                ra.set_additional_info("FSM_class", "C")
                self.transcript2transcript.append(ra)
                continue

            fsm_class = "A" if len(gene_to_model[ra.isoform_matches[0].assigned_gene]) == 1 else "C"
            ra.set_additional_info("FSM_class", fsm_class)
            assigned_tid = ra.isoform_matches[0].assigned_transcript
            if not assigned_tid or assigned_tid not in self.gene_info.all_isoforms_introns:
                continue

            ref_introns = self.gene_info.all_isoforms_introns[assigned_tid]
            iso_introns = junctions_from_blocks(model.exon_blocks)
            event_str = ",".join(
                match_subtype_to_str_with_additional_info(x, model.strand, iso_introns, ref_introns)
                for x in ra.isoform_matches[0].match_subclassifications
            )
            model.add_additional_attribute("similar_reference_id", assigned_tid)
            model.add_additional_attribute("alternatives", event_str)
            self.transcript2transcript.append(ra)
