"""Stage 1 of transcript discovery: full-length isoforms from IntronGraph paths."""

from __future__ import annotations

import logging
from functools import cmp_to_key
from typing import Optional, Tuple

from isoquant_next.alignment.alignment_info import PolyAInfo
from isoquant_next.annotation.gene_info import TranscriptModel, TranscriptModelType
from isoquant_next.assignment.types import ReadAssignment, is_matching_assignment
from isoquant_next.core.utils import TranscriptNaming, cmp, get_exons
from isoquant_next.discovery.intron_graph import TerminalVertex
from isoquant_next.discovery.state import StrandnessReportingLevel

logger = logging.getLogger("IsoQuant")


class GraphBasedConstructor:
    """Builds full-length transcript models from IntronGraph FL paths."""

    def __init__(self, ctx) -> None:
        self.ctx = ctx
        self.store = ctx.store
        self.args = ctx.args
        self.gene_info = ctx.gene_info
        self.use_tss_model = ctx.use_tss_model
        self.create_alt_ts_nics = ctx.create_alt_ts_nics
        self.assigner = ctx.assigner
        self.profile_constructor = ctx.profile_constructor
        self.strand_detector = ctx.strand_detector
        self.use_technical_replicas = ctx.use_technical_replicas
        self.file_name_group_idx = ctx.file_name_group_idx

    def _reference_isoform_for_path(
        self, assignment: ReadAssignment, path: tuple, intron_path: tuple, transcript_range: Tuple[int, int]
    ) -> Optional[str]:
        if is_matching_assignment(assignment):
            matched_ref = assignment.isoform_matches[0].assigned_transcript
        elif intron_path in self.ctx.known_isoforms_in_graph:
            matched_ref = self.ctx.known_isoforms_in_graph[intron_path]
        else:
            return None

        if matched_ref is None or matched_ref not in self.gene_info.all_isoforms_exons:
            return None

        ref_exons = self.gene_info.all_isoforms_exons[matched_ref]
        left_diff = abs(transcript_range[0] - ref_exons[0][0]) > self.args.apa_delta
        right_diff = abs(transcript_range[1] - ref_exons[-1][1]) > self.args.apa_delta

        if path[0][0] == TerminalVertex.polyt and left_diff:
            return None
        if path[-1][0] == TerminalVertex.polya and right_diff:
            return None

        if self.use_tss_model:
            strand = self.gene_info.isoform_strands.get(matched_ref, ".")
            if strand == "+" and path[0][0] == TerminalVertex.tss_left and left_diff:
                return None
            if strand == "-" and path[-1][0] == TerminalVertex.tss_right and right_diff:
                return None

        return matched_ref

    def construct_fl_isoforms(self) -> None:
        logger.debug("Total FL paths %d", len(self.ctx.path_storage.fl_paths))
        ordered_paths = sorted(
            self.ctx.path_storage.fl_paths,
            key=cmp_to_key(lambda x, y: cmp(x, y) if len(x) == len(y) else cmp(len(y), len(x))),
        )
        for path in ordered_paths:
            intron_path = path[1:-1]
            if not intron_path:
                continue
            transcript_range = (path[0][1], path[-1][1])
            novel_exons = get_exons(transcript_range, list(intron_path))
            count = self.ctx.path_storage.paths[path]
            new_tid = TranscriptNaming.transcript_prefix + str(self.ctx.get_transcript_id())

            reference_isoform = None
            if intron_path[0][0] == TerminalVertex.polyt:
                polya_info = PolyAInfo(-1, intron_path[0][1], -1, -1)
            elif intron_path[-1][0] == TerminalVertex.polya:
                polya_info = PolyAInfo(intron_path[-1][1], -1, -1, -1)
            else:
                polya_info = PolyAInfo(-1, -1, -1, -1)

            combined_profile = self.profile_constructor.construct_profiles(novel_exons, polya_info, [])
            assignment = self.assigner.assign_to_isoform(new_tid, combined_profile)

            if self.create_alt_ts_nics:
                reference_isoform = self._reference_isoform_for_path(assignment, path, intron_path, transcript_range)
            elif is_matching_assignment(assignment):
                reference_isoform = assignment.isoform_matches[0].assigned_transcript
            elif intron_path in self.ctx.known_isoforms_in_graph:
                continue

            new_model = None
            if reference_isoform:
                if reference_isoform not in self.store.detected_known_isoforms and count >= self.args.min_known_count:
                    new_model = self.store.transcript_from_reference(reference_isoform)
                    self.store.detected_known_isoforms.add(reference_isoform)
            else:
                has_polyt = path[0][0] == TerminalVertex.polyt
                has_polya = path[-1][0] == TerminalVertex.polya
                polya_site = has_polya or has_polyt
                tx_strand = self.strand_detector.get_strand(intron_path, has_polya, has_polyt)
                clean_strand = self.strand_detector.get_clean_strand(intron_path)

                if count < self.args.min_novel_count:
                    pass
                elif len(novel_exons) == 2 and (
                    (self.args.require_monointronic_polya and not polya_site) or clean_strand == "."
                ):
                    pass
                elif (
                    self.args.report_canonical_strategy == StrandnessReportingLevel.only_canonical and clean_strand == "."
                ) or (
                    self.args.report_canonical_strategy == StrandnessReportingLevel.only_stranded and tx_strand == "."
                ):
                    logger.debug("Avoiding unreliable transcript with %d exons (strand unknown)", len(novel_exons))
                else:
                    if self.use_technical_replicas and self.file_name_group_idx >= 0:
                        read_assignments = self.ctx.path_storage.paths_to_reads[path]
                        if read_assignments:
                            file_groups = {
                                a.read_group[self.file_name_group_idx]
                                for a in read_assignments
                                if self.file_name_group_idx < len(a.read_group)
                            }
                            if len(file_groups) <= 1:
                                continue

                    tx_gene = self.ctx.select_reference_gene(intron_path, transcript_range, tx_strand)
                    if tx_gene is None:
                        tx_gene = f"{TranscriptNaming.novel_gene_prefix}{self.gene_info.chr_id}_{self.ctx.get_transcript_id()}"
                    elif tx_strand == ".":
                        tx_strand = self.gene_info.gene_strands[tx_gene]

                    if all(intron in self.ctx.known_introns for intron in intron_path):
                        tx_type = TranscriptModelType.novel_in_catalog
                        id_suffix = TranscriptNaming.nic_transcript_suffix
                    else:
                        tx_type = TranscriptModelType.novel_not_in_catalog
                        id_suffix = TranscriptNaming.nnic_transcript_suffix

                    new_model = TranscriptModel(
                        self.gene_info.chr_id,
                        tx_strand,
                        f"{new_tid}.{self.gene_info.chr_id}{id_suffix}",
                        tx_gene,
                        novel_exons,
                        tx_type,
                    )
                    new_model.intron_path = intron_path

            if new_model is not None:
                self.store.transcript_model_storage.append(new_model)
                for ra in self.ctx.path_storage.paths_to_reads[path]:
                    self.store.save_assigned_read(ra, new_model.transcript_id)
