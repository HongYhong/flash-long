"""Fusion candidate metadata construction and breakpoint consensus."""

from __future__ import annotations

import logging
from typing import List, Optional

logger = logging.getLogger("IsoQuant")


class FusionMetadata:
    """Build per-fusion consensus metadata, gene assignments, and biotypes."""

    def __init__(self, detector) -> None:
        self.detector = detector

    def process_all(self, min_support: int = 1) -> None:
        detector = self.detector
        logger.info(
            "FusionMetadata: Processing %d fusion keys with min_support=%d",
            len(detector.fusion_candidates),
            min_support,
        )
        for fusion_key, reads in list(detector.fusion_candidates.items()):
            self._process_fusion_candidate(fusion_key, reads, min_support)

    def _process_fusion_candidate(self, fusion_key: str, reads: set, min_support: int) -> None:
        detector = self.detector
        gate_result = self._passes_all_gates(fusion_key, reads, min_support)
        if not gate_result:
            meta = detector.fusion_metadata.get(fusion_key)
            if meta is not None:
                meta["is_valid"] = False
                meta.setdefault("reasons", []).append("Failed gating (support or clustering)")
            return
        _support, consensus_bp, clustered_support = gate_result
        left_chr, left_pos, right_chr, right_pos = consensus_bp
        meta = self._initialize_or_get_metadata(fusion_key, reads, clustered_support, consensus_bp)
        self._assign_genes(fusion_key, meta, left_chr, left_pos, right_chr, right_pos)
        left_gene = meta.get("left_gene")
        right_gene = meta.get("right_gene")
        self._update_fusion_key_mappings(fusion_key, left_gene, right_gene, meta)
        self._compute_biotypes(meta, left_gene, right_gene, left_chr, left_pos, right_chr, right_pos)
        meta["is_valid"] = True
        meta.setdefault("reasons", [])

    def _initialize_or_get_metadata(
        self,
        fusion_key: str,
        reads: set,
        clustered_support: int,
        consensus_bp: tuple,
    ) -> dict:
        detector = self.detector
        meta = detector.fusion_metadata.get(fusion_key)
        if meta is None:
            meta = {"supporting_reads": set(), "left_gene": None, "right_gene": None}
            detector.fusion_metadata[fusion_key] = meta
        meta["supporting_reads"] = set(reads)
        meta["support"] = clustered_support
        meta["consensus_bp"] = consensus_bp
        return meta

    def _compute_biotypes(
        self,
        meta: dict,
        left_gene: Optional[str],
        right_gene: Optional[str],
        left_chr: str,
        left_pos: int,
        right_chr: str,
        right_pos: int,
    ) -> None:
        detector = self.detector
        if left_gene:
            meta["left_biotype"] = detector.get_gene_biotype(left_gene, chrom=left_chr, pos=left_pos)
        if right_gene:
            meta["right_biotype"] = detector.get_gene_biotype(right_gene, chrom=right_chr, pos=right_pos)

    def _update_fusion_key_mappings(
        self,
        original_key: str,
        final_left_gene: Optional[str],
        final_right_gene: Optional[str],
        meta: dict,
    ) -> str:
        detector = self.detector
        new_left = final_left_gene or "intergenic"
        new_right = final_right_gene or "intergenic"
        sorted_genes = sorted([new_left, new_right])
        new_key = f"{sorted_genes[0]}--{sorted_genes[1]}"
        if new_key != original_key:
            logger.info("Gene assignment changed: %s -> %s", original_key, new_key)
            detector.fusion_metadata[new_key] = meta
            detector.fusion_metadata.pop(original_key, None)
            detector.fusion_candidates[new_key] = detector.fusion_candidates.pop(original_key, set())
            if original_key in detector.fusion_breakpoints:
                detector.fusion_breakpoints[new_key] = detector.fusion_breakpoints.pop(original_key)
            if original_key in detector.fusion_assigned_pairs:
                detector.fusion_assigned_pairs[new_key] = detector.fusion_assigned_pairs.pop(original_key)
        return new_key

    def _passes_all_gates(self, fusion_key: str, reads: set, min_support: int):
        support = len(reads)
        if support < min_support:
            return None
        bp_counts = self.detector.fusion_breakpoints.get(fusion_key, {})
        if not bp_counts:
            return None
        consensus_result = self.detector.cluster_breakpoints(bp_counts, window=2000)
        if not consensus_result:
            return None
        consensus_bp, clustered_support = consensus_result
        return support, consensus_bp, clustered_support

    def _assign_genes(
        self,
        fusion_key: str,
        meta: dict,
        left_chr: str,
        left_pos: int,
        right_chr: str,
        right_pos: int,
    ) -> None:
        parts = fusion_key.split("--")
        if len(parts) == 2:
            left_gene, right_gene = parts[0], parts[1]
        else:
            left_gene, right_gene = None, None
        meta["left_gene"] = left_gene
        meta["right_gene"] = right_gene
        read_scores = self.detector.fusion_read_scores.get(fusion_key, {})
        left_scores: List[float] = []
        right_scores: List[float] = []
        for read_name in meta["supporting_reads"]:
            if read_name in read_scores:
                ls, rs = read_scores[read_name]
                left_scores.append(ls)
                right_scores.append(rs)
        meta["left_score"] = sum(left_scores) / len(left_scores) if left_scores else 0.0
        meta["right_score"] = sum(right_scores) / len(right_scores) if right_scores else 0.0
