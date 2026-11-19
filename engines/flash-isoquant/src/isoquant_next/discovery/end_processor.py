"""Stage 3 of transcript discovery: terminal end refinement & alternative-end NIC discovery."""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Callable, Dict, List, Optional, Tuple

from isoquant_next.annotation.gene_info import TranscriptModel, TranscriptModelType
from isoquant_next.assignment.types import ReadAssignment
from isoquant_next.core.utils import TranscriptNaming, junctions_from_blocks
from isoquant_next.terminal.peaks import detect_peaks, get_polya_model, get_tss_model

logger = logging.getLogger("IsoQuant")


class TranscriptEndProcessor:
    """Refines transcript terminal coordinates via L1 GBDT peaks & creates alternative-end NICs."""

    def __init__(self, ctx) -> None:
        self.ctx = ctx
        self.store = ctx.store
        self.gene_info = ctx.gene_info
        self.args = ctx.args
        self.use_tss_model = ctx.use_tss_model
        self.novel_apa = ctx.novel_apa

    def correct_novel_transcript_ends(
        self, transcript_model: TranscriptModel, assigned_reads: List[ReadAssignment]
    ) -> None:
        transcript_start = transcript_model.exon_blocks[0][0]
        transcript_end = transcript_model.exon_blocks[-1][1]
        first_exon_right = transcript_model.exon_blocks[0][1]
        last_exon_left = transcript_model.exon_blocks[-1][0]
        strand = transcript_model.strand

        start_reads, end_reads = self._terminal_histograms(transcript_model, assigned_reads)
        start_hist = {p: len(rs) for p, rs in start_reads.items()}
        end_hist = {p: len(rs) for p, rs in end_reads.items()}
        start_supported = any(abs(p - transcript_start) <= self.args.apa_delta for p in start_hist)
        end_supported = any(abs(p - transcript_end) <= self.args.apa_delta for p in end_hist)

        if not start_supported:
            model = self._terminal_model(strand, left=True)
            if model is not None:
                new_start = self._peak_boundary(start_hist, model, lambda pos: pos < first_exon_right)
            else:
                new_start = self._closest_inward(sorted(start_hist.keys()), transcript_start, greater=True)
            if new_start is not None and new_start != transcript_start and new_start < first_exon_right:
                transcript_model.exon_blocks[0] = (new_start, first_exon_right)

        if not end_supported:
            model = self._terminal_model(strand, left=False)
            if model is not None:
                new_end = self._peak_boundary(end_hist, model, lambda pos: pos > last_exon_left)
            else:
                new_end = self._closest_inward(sorted(end_hist.keys(), reverse=True), transcript_end, greater=False)
            if new_end is not None and new_end != transcript_end and new_end > last_exon_left:
                transcript_model.exon_blocks[-1] = (last_exon_left, new_end)

    def _mark_terminal_confirmation(self, model: TranscriptModel) -> None:
        reads = self.store.transcript_read_ids[model.transcript_id]
        if not reads:
            model.polya_confirmed = model.tss_confirmed = False
            return
        start_reads, end_reads = self._terminal_histograms(model, reads)
        m_start = model.exon_blocks[0][0]
        m_end = model.exon_blocks[-1][1]
        start_ok = any(abs(p - m_start) <= self.args.apa_delta for p in start_reads)
        end_ok = any(abs(p - m_end) <= self.args.apa_delta for p in end_reads)
        if model.strand == "+":
            model.tss_confirmed, model.polya_confirmed = start_ok, end_ok
        elif model.strand == "-":
            model.polya_confirmed, model.tss_confirmed = start_ok, end_ok
        else:
            model.polya_confirmed = model.tss_confirmed = False

    def _terminal_model(self, strand: str, left: bool):
        if strand == "+":
            is_polya_side = not left
        elif strand == "-":
            is_polya_side = left
        else:
            return None
        if is_polya_side:
            return get_polya_model()
        return get_tss_model() if self.use_tss_model else None

    @staticmethod
    def _peak_boundary(histogram: Dict[int, int], model, valid: Callable[[int], bool]) -> Optional[int]:
        if not histogram:
            return None
        peaks = [p for p in detect_peaks(histogram, model) if valid(p.position)]
        if not peaks:
            return None
        return max(peaks, key=lambda p: p.count).position

    @staticmethod
    def _closest_inward(sorted_positions: List[int], boundary: int, greater: bool) -> Optional[int]:
        for pos in sorted_positions:
            if greater and pos > boundary:
                return pos
            if not greater and pos < boundary:
                return pos
        return None

    @staticmethod
    def _confirmed_polya_pos(assignment: ReadAssignment, strand: str) -> Optional[int]:
        info = getattr(assignment, "polya_info", None)
        if not getattr(assignment, "polyA_found", False) or info is None:
            return None
        if strand == "+" and info.external_polya_pos != -1:
            return info.external_polya_pos
        if strand == "-" and info.external_polyt_pos != -1:
            return info.external_polyt_pos
        return None

    def _terminal_histograms(
        self, source_model: TranscriptModel, assigned_reads: List[ReadAssignment]
    ) -> Tuple[Dict[int, List[ReadAssignment]], Dict[int, List[ReadAssignment]]]:
        strand = source_model.strand
        first_right = source_model.exon_blocks[0][1]
        last_left = source_model.exon_blocks[-1][0]
        start_reads: Dict[int, List[ReadAssignment]] = defaultdict(list)
        end_reads: Dict[int, List[ReadAssignment]] = defaultdict(list)

        for a in assigned_reads:
            ex = a.corrected_exons
            if strand == "+":
                if ex[0][0] < first_right:
                    start_reads[ex[0][0]].append(a)
                pos = self._confirmed_polya_pos(a, "+")
                if pos is not None and pos > last_left:
                    end_reads[pos].append(a)
            elif strand == "-":
                pos = self._confirmed_polya_pos(a, "-")
                if pos is not None and pos < first_right:
                    start_reads[pos].append(a)
                if ex[-1][1] > last_left:
                    end_reads[ex[-1][1]].append(a)
            else:
                if ex[0][0] < first_right:
                    start_reads[ex[0][0]].append(a)
                if ex[-1][1] > last_left:
                    end_reads[ex[-1][1]].append(a)
        return start_reads, end_reads

    @staticmethod
    def _intron_chain_key(model: TranscriptModel) -> Tuple[Tuple[int, int], ...]:
        return tuple(junctions_from_blocks(model.exon_blocks))

    def _add_known_alternative_end_models(self) -> None:
        existing_pairs = defaultdict(list)
        for ref_exons in self.gene_info.all_isoforms_exons.values():
            ck = tuple(junctions_from_blocks(ref_exons))
            existing_pairs[ck].append((ref_exons[0][0], ref_exons[-1][1]))
        for m in self.store.transcript_model_storage:
            existing_pairs[self._intron_chain_key(m)].append((m.exon_blocks[0][0], m.exon_blocks[-1][1]))

        new_models = []
        for model in self.store.transcript_model_storage:
            if model.transcript_type != TranscriptModelType.known and not self.novel_apa:
                continue
            source_id = model.transcript_id
            for nic, nic_reads in self.derive_alternative_end_models(model, self.store.transcript_read_ids[source_id]):
                ck = self._intron_chain_key(nic)
                ns, ne = nic.exon_blocks[0][0], nic.exon_blocks[-1][1]
                if any(
                    abs(ns - s) <= self.args.apa_delta and abs(ne - e) <= self.args.apa_delta
                    for s, e in existing_pairs[ck]
                ):
                    continue
                existing_pairs[ck].append((ns, ne))
                new_models.append(nic)
                for read in nic_reads:
                    if len(self.store.transcript_read_ids[source_id]) <= 1:
                        break
                    self.store.move_read_to_model(read, source_id, nic.transcript_id)
        if new_models:
            self.store.transcript_model_storage.extend(new_models)

    def _drop_duplicate_alt_end_models(self) -> None:
        ref_pairs = defaultdict(list)
        for ref_exons in self.gene_info.all_isoforms_exons.values():
            ref_pairs[tuple(junctions_from_blocks(ref_exons))].append((ref_exons[0][0], ref_exons[-1][1]))

        d = self.args.apa_delta
        kept = []
        kept_pairs = defaultdict(list)
        for m in self.store.transcript_model_storage:
            if m.transcript_type == TranscriptModelType.known:
                kept.append(m)
                continue
            ck = self._intron_chain_key(m)
            s, e = m.exon_blocks[0][0], m.exon_blocks[-1][1]
            if any(abs(s - rs) <= d and abs(e - re) <= d for rs, re in ref_pairs[ck]) or any(
                abs(s - ks) <= d and abs(e - ke) <= d for ks, ke in kept_pairs[ck]
            ):
                if m.transcript_id in self.store.internal_counter:
                    self.store.delete_from_storage(m.transcript_id)
                continue
            kept.append(m)
            kept_pairs[ck].append((s, e))
        self.store.transcript_model_storage = kept

    def derive_alternative_end_models(
        self, source_model: TranscriptModel, assigned_reads: List[ReadAssignment]
    ) -> List[Tuple[TranscriptModel, List[ReadAssignment]]]:
        if not assigned_reads:
            return []
        exon_blocks = source_model.exon_blocks
        strand = source_model.strand
        ann_start = exon_blocks[0][0]
        ann_end = exon_blocks[-1][1]
        first_right = exon_blocks[0][1]
        last_left = exon_blocks[-1][0]
        sibling_type = (
            TranscriptModelType.novel_in_catalog
            if source_model.transcript_type == TranscriptModelType.known
            else source_model.transcript_type
        )

        start_reads, end_reads = self._terminal_histograms(source_model, assigned_reads)
        start_hist = {p: len(rs) for p, rs in start_reads.items()}
        end_hist = {p: len(rs) for p, rs in end_reads.items()}
        claimed: set[int] = set()

        def assign_to_peaks(reads_map, alt_peaks, annotated_pos):
            per_peak = {p: [] for p in alt_peaks}
            candidates = list(alt_peaks) + [annotated_pos]
            for p, rs in reads_map.items():
                nearest = min(candidates, key=lambda c: (abs(c - p), c))
                if nearest != annotated_pos and abs(nearest - p) <= self.args.apa_delta:
                    for r in rs:
                        if id(r) not in claimed:
                            claimed.add(id(r))
                            per_peak[nearest].append(r)
            return per_peak

        new_models: List[Tuple[TranscriptModel, List[ReadAssignment]]] = []
        start_model = self._terminal_model(strand, left=True)
        if start_model is not None:
            start_peaks = self._alternative_end_positions(
                start_hist, start_model, ann_start, lambda p: p < first_right
            )
            per_peak = assign_to_peaks(start_reads, start_peaks, ann_start)
            for pos in start_peaks:
                new_models.append((self._nic_model_with_boundary(source_model, sibling_type, start=pos), per_peak[pos]))

        end_model = self._terminal_model(strand, left=False)
        if end_model is not None:
            end_peaks = self._alternative_end_positions(
                end_hist, end_model, ann_end, lambda p: p > last_left
            )
            per_peak = assign_to_peaks(end_reads, end_peaks, ann_end)
            for pos in end_peaks:
                new_models.append((self._nic_model_with_boundary(source_model, sibling_type, end=pos), per_peak[pos]))
        return new_models

    def _alternative_end_positions(
        self, histogram: Dict[int, int], model, annotated_pos: int, clamp: Callable[[int], bool]
    ) -> List[int]:
        if not histogram:
            return []
        peaks = [p for p in detect_peaks(histogram, model) if clamp(p.position)]
        if not peaks:
            return []
        dominant = max(p.count for p in peaks)
        cutoff = max(self.args.min_novel_count, self.args.terminal_position_rel * dominant)
        return [
            p.position
            for p in peaks
            if abs(p.position - annotated_pos) > self.args.apa_delta and p.count >= cutoff
        ]

    def _nic_model_with_boundary(
        self,
        source_model: TranscriptModel,
        transcript_type: TranscriptModelType,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> TranscriptModel:
        exon_blocks = [tuple(e) for e in source_model.exon_blocks]
        if start is not None:
            exon_blocks[0] = (start, exon_blocks[0][1])
        if end is not None:
            exon_blocks[-1] = (exon_blocks[-1][0], end)
        suffix = (
            TranscriptNaming.nic_transcript_suffix
            if transcript_type == TranscriptModelType.novel_in_catalog
            else TranscriptNaming.nnic_transcript_suffix
        )
        new_tid = TranscriptNaming.transcript_prefix + str(self.ctx.get_transcript_id())
        new_model = TranscriptModel(
            self.gene_info.chr_id,
            source_model.strand,
            f"{new_tid}.{self.gene_info.chr_id}{suffix}",
            source_model.gene_id,
            exon_blocks,
            transcript_type,
        )
        new_model.intron_path = source_model.intron_path
        return new_model
