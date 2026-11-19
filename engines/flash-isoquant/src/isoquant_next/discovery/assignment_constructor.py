"""Stage 2 of transcript discovery: non-FL known, mono-exon known, and novel mono-exon models."""

from __future__ import annotations

import logging
from collections import defaultdict

from isoquant_next.annotation.gene_info import TranscriptModel, TranscriptModelType
from isoquant_next.assignment.types import MatchEventSubtype
from isoquant_next.core.utils import TranscriptNaming, intersection_len, interval_len
from isoquant_next.discovery.intron_graph import TerminalVertex

logger = logging.getLogger("IsoQuant")


class AssignmentBasedConstructor:
    """Builds known non-FL and mono-exon models from per-read reference assignments."""

    def __init__(self, ctx) -> None:
        self.ctx = ctx
        self.store = ctx.store
        self.args = ctx.args
        self.gene_info = ctx.gene_info

    def construct_assignment_based_isoforms(self, read_assignment_storage) -> None:
        spliced_reads = defaultdict(list)
        left_support = defaultdict(int)
        right_support = defaultdict(int)
        polya_sites = defaultdict(int)
        mono_reads = defaultdict(list)
        mono_cov: dict[str, bytearray] = {}
        novel_mono_reads = []

        for ra in read_assignment_storage:
            if not ra:
                continue
            if len(ra.corrected_exons) <= 2 and (
                ra.multimapper or ra.mapping_quality < self.args.simple_alignments_mapq_cutoff
            ):
                continue

            if (
                len(ra.corrected_exons) == 1
                and ra.polyA_found
                and not ra.multimapper
                and (ra.assignment_type.is_unassigned() or ra.assignment_type.is_inconsistent())
            ):
                novel_mono_reads.append(ra)

            if not ra.assignment_type.is_unique():
                continue

            ref_id = ra.isoform_matches[0].assigned_transcript
            if ref_id in self.store.detected_known_isoforms:
                continue

            events = ra.isoform_matches[0].match_subclassifications
            if any(e.event_type == MatchEventSubtype.mono_exon_match for e in events):
                mono_reads[ref_id].append(ra)
                assert len(self.gene_info.all_isoforms_introns[ref_id]) == 0
                t_start = self.gene_info.all_isoforms_exons[ref_id][0][0]
                t_end = self.gene_info.all_isoforms_exons[ref_id][-1][1]
                t_len = t_end - t_start + 1
                cov = mono_cov.get(ref_id)
                if cov is None:
                    cov = bytearray(t_len)
                    mono_cov[ref_id] = cov
                s = max(0, ra.corrected_exons[0][0] - t_start)
                e = min(t_len, ra.corrected_exons[-1][1] - t_start + 1)
                if e > s:
                    cov[s:e] = b"\x01" * (e - s)

                if self.gene_info.isoform_strands[ref_id] == "+":
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_right for x in events):
                        polya_sites[ref_id] += 1
                else:
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_left for x in events):
                        polya_sites[ref_id] += 1
            elif len(self.gene_info.all_isoforms_exons[ref_id]) > 1:
                spliced_reads[ref_id].append(ra)
                if self.args.requires_polya_for_construction and self.gene_info.isoform_strands[ref_id] == "-":
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_left for x in events):
                        left_support[ref_id] += 1
                elif abs(self.gene_info.all_isoforms_exons[ref_id][0][0] - ra.corrected_exons[0][0]) <= self.args.apa_delta:
                    left_support[ref_id] += 1

                if self.args.requires_polya_for_construction and self.gene_info.isoform_strands[ref_id] == "+":
                    if any(x.event_type == MatchEventSubtype.correct_polya_site_right for x in events):
                        right_support[ref_id] += 1
                elif abs(self.gene_info.all_isoforms_exons[ref_id][-1][1] - ra.corrected_exons[-1][1]) <= self.args.apa_delta:
                    right_support[ref_id] += 1

        self.construct_monoexon_isoforms(mono_reads, mono_cov, polya_sites)
        if not self.args.fl_only:
            self.construct_nonfl_isoforms(spliced_reads, left_support, right_support)
        if self.args.report_novel_unspliced:
            self.construct_monoexon_novel(novel_mono_reads)

    def collect_terminal_exons_from_graph(self):
        polya_exons = []
        polyt_exons = []
        for intron, targets in self.ctx.intron_graph.outgoing_edges.items():
            for v in targets:
                if v[0] == TerminalVertex.polya:
                    polya_exons.append((intron[1], v[1]))
        for intron, sources in self.ctx.intron_graph.incoming_edges.items():
            for v in sources:
                if v[0] == TerminalVertex.polyt:
                    polyt_exons.append((v[1], intron[0]))
        return polya_exons, polyt_exons

    def is_internal_monoexonic_read(self, alignment, terminal_exons, forward: bool = True) -> bool:
        read_coord = alignment.corrected_exons[0]
        if forward:
            return any(
                abs(e[1] - alignment.corrected_exons[-1][1]) <= self.args.apa_delta
                and read_coord[0] >= e[0] - self.args.delta
                for e in terminal_exons
            )
        return any(
            abs(e[0] - alignment.corrected_exons[0][0]) <= self.args.apa_delta
            and read_coord[1] <= e[1] + self.args.delta
            for e in terminal_exons
        )

    def construct_monoexon_novel(self, novel_mono_reads) -> None:
        polya_exons, polyt_exons = self.collect_terminal_exons_from_graph()
        polya_reads = defaultdict(list)
        polyt_reads = defaultdict(list)
        for a in novel_mono_reads:
            if a.polya_info.external_polya_pos != -1 and not self.is_internal_monoexonic_read(a, polya_exons, forward=True):
                polya_reads[a.polya_info.external_polya_pos].append(a)
            if a.polya_info.external_polyt_pos != -1 and not self.is_internal_monoexonic_read(a, polyt_exons, forward=False):
                polyt_reads[a.polya_info.external_polyt_pos].append(a)

        self.generate_monoexon_from_clustered(self.cluster_monoexons(polya_reads), True)
        self.generate_monoexon_from_clustered(self.cluster_monoexons(polyt_reads), False)

    def generate_monoexon_from_clustered(self, clustered_reads, forward: bool = True) -> set:
        cutoff = self.args.min_novel_count
        result = set()
        for three_prime_pos, reads in clustered_reads.items():
            if len(reads) < cutoff:
                continue
            if forward:
                five_prime_pos = min(a.corrected_exons[0][0] for a in reads)
                coordinates = (five_prime_pos, three_prime_pos)
                strand = "+"
            else:
                five_prime_pos = max(a.corrected_exons[-1][1] for a in reads)
                coordinates = (three_prime_pos, five_prime_pos)
                strand = "-"

            new_tid = TranscriptNaming.transcript_prefix + str(self.ctx.get_transcript_id())
            tx_gene = f"{TranscriptNaming.novel_gene_prefix}{self.gene_info.chr_id}_{self.ctx.get_transcript_id()}"
            half_len = interval_len(coordinates) / 2
            if any(
                any(intersection_len(exon, coordinates) > half_len for exon in m.exon_blocks)
                for m in self.store.transcript_model_storage
            ):
                continue

            new_model = TranscriptModel(
                self.gene_info.chr_id,
                strand,
                f"{new_tid}.{self.gene_info.chr_id}{TranscriptNaming.nnic_transcript_suffix}",
                tx_gene,
                [coordinates],
                TranscriptModelType.novel_not_in_catalog,
            )
            result.add(coordinates)
            self.store.transcript_model_storage.append(new_model)
            for ra in reads:
                self.store.save_assigned_read(ra, new_model.transcript_id)
        return result

    def cluster_monoexons(self, grouped_reads):
        clustered = defaultdict(list)
        while grouped_reads:
            top_pos = max(grouped_reads.items(), key=lambda x: len(x[1]))[0]
            for pos in range(top_pos - self.args.apa_delta, top_pos + self.args.apa_delta + 1):
                if pos in grouped_reads:
                    clustered[top_pos] += grouped_reads.pop(pos)
        return clustered

    def construct_monoexon_isoforms(self, mono_reads, mono_cov, polya_sites) -> None:
        for isoform_id, reads in mono_reads.items():
            cov = mono_cov[isoform_id]
            coverage = float(cov.count(1)) / float(len(cov))
            polya_support = polya_sites[isoform_id]
            if (
                len(reads) < self.args.min_known_count
                or coverage < self.args.min_mono_exon_coverage
                or (self.args.require_monoexonic_polya and polya_support == 0)
            ):
                continue
            if isoform_id not in self.store.detected_known_isoforms:
                new_model = self.store.transcript_from_reference(isoform_id)
                self.store.transcript_model_storage.append(new_model)
                self.store.detected_known_isoforms.add(isoform_id)
                for ra in reads:
                    self.store.save_assigned_read(ra, new_model.transcript_id)

    def construct_nonfl_isoforms(self, spliced_reads, left_support, right_support) -> None:
        for isoform_id, reads in spliced_reads.items():
            if isoform_id in self.store.detected_known_isoforms:
                continue
            if isoform_id not in self.ctx.known_isoforms_in_graph_ids:
                continue
            if len(reads) < self.args.min_known_count or left_support[isoform_id] < 1 or right_support[isoform_id] < 1:
                continue
            new_model = self.store.transcript_from_reference(isoform_id)
            self.store.transcript_model_storage.append(new_model)
            self.store.detected_known_isoforms.add(isoform_id)
            self.store.nonfl_known_reads[new_model.transcript_id] = list(reads)
            for ra in reads:
                self.store.save_assigned_read(ra, new_model.transcript_id)
