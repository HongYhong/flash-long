"""Fusion candidate filtering, duplicate merging, and confidence scoring."""

from __future__ import annotations

from collections import defaultdict
import logging
from typing import Optional

logger = logging.getLogger("IsoQuant")

ALLOWED_BIOTYPES = frozenset({
    "protein_coding",
    "IG_C_gene",
    "IG_V_gene",
    "IG_D_gene",
    "IG_J_gene",
    "TR_V_gene",
    "TR_D_gene",
    "TR_J_gene",
    "TR_C_gene",
})

DRIVER_GENES = frozenset({
    "BCR", "ABL1", "RUNX1", "RUNX1T1", "PML", "RARA", "ETV6",
    "NUP98", "NUP214", "KMT2A", "MLLT3", "EWSR1", "ERG",
    "FGFR1", "RET", "ROS1", "ALK", "TMPRSS2", "FLI1",
    "MYH11", "CBFB",
})

MULTICOPY_PREFIXES = (
    "RPL", "RPS", "MRPL", "MRPS", "H1-", "H2A", "H2B", "H3-", "H4-",
    "HIST1", "HIST2", "HIST3", "HLA-", "MICA", "MICB", "KRT",
    "OR", "ZNF", "DEFA", "DEFB", "IGH", "IGK", "IGL", "TRAV", "TRBV", "TRGV", "TRDV",
    "MUC", "AMY", "CYP",
)

RIBOSOMAL_HISTONE_PREFIXES = ("RPL", "RPS", "RPLP", "HIST", "H1", "H2A", "H2B", "H3", "H4")


class FusionValidator:
    """Stateless validator operating over a FusionDetector's candidate state."""

    def __init__(self, detector) -> None:
        self.detector = detector

    def _get_breakpoint_coords(self, fusion_key: str):
        left_chr = left_pos = right_chr = right_pos = None
        bp_counts = self.detector.fusion_breakpoints.get(fusion_key, {})
        if bp_counts:
            first_bp = next(iter(bp_counts.keys()), None)
            if first_bp and len(first_bp) == 4:
                left_chr, left_pos, right_chr, right_pos = first_bp
        return left_chr, left_pos, right_chr, right_pos

    def _is_allowed_biotype(self, biotype: Optional[str]) -> bool:
        return biotype in ALLOWED_BIOTYPES

    def _salvage_and_check_gene_pair(
        self, left_gene, right_gene, left_chr, left_pos, right_chr, right_pos
    ):
        left_biotype = self.detector.get_gene_biotype(left_gene, chrom=left_chr, pos=left_pos)
        right_biotype = self.detector.get_gene_biotype(right_gene, chrom=right_chr, pos=right_pos)
        has_non_coding = False
        reason = None
        if left_gene and left_biotype and not self._is_allowed_biotype(left_biotype):
            has_non_coding = True
            reason = f"{left_gene}={left_biotype}"
        elif right_gene and right_biotype and not self._is_allowed_biotype(right_biotype):
            has_non_coding = True
            reason = f"{right_gene}={right_biotype}"
        return left_gene, left_biotype, right_gene, right_biotype, has_non_coding, reason

    def filter_early_non_coding_genes(self) -> None:
        fusions_to_discard = set()
        for fusion_key, read_assignments in self.detector.fusion_assigned_pairs.items():
            if fusion_key not in self.detector.fusion_candidates:
                continue
            left_chr, left_pos, right_chr, right_pos = self._get_breakpoint_coords(fusion_key)
            has_non_coding = False
            for read_name, (left_gene, right_gene) in read_assignments.items():
                left_gene, _lb, right_gene, _rb, is_nc, _reason = self._salvage_and_check_gene_pair(
                    left_gene, right_gene, left_chr, left_pos, right_chr, right_pos
                )
                if self.detector.fusion_assigned_pairs.get(fusion_key, {}).get(read_name):
                    self.detector.fusion_assigned_pairs[fusion_key][read_name] = (left_gene, right_gene)
                if is_nc:
                    has_non_coding = True
                    break
            if has_non_coding:
                fusions_to_discard.add(fusion_key)
        self._remove_discarded_fusions_internal(fusions_to_discard)

    def _remove_discarded_fusions_internal(self, fusions_to_discard) -> None:
        for fusion_key in fusions_to_discard:
            self.detector.fusion_metadata.pop(fusion_key, None)
            self.detector.fusion_candidates.pop(fusion_key, None)
            self.detector.fusion_breakpoints.pop(fusion_key, None)
            self.detector.fusion_assigned_pairs.pop(fusion_key, None)

    def filter_based_on_biotype(self) -> None:
        fusions_to_discard = []
        for fusion_key, meta in list(self.detector.fusion_metadata.items()):
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            left_chr = left_pos = right_chr = right_pos = None
            consensus_bp = meta.get("consensus_bp")
            if consensus_bp:
                left_chr, left_pos, right_chr, right_pos = consensus_bp
            left_biotype = self.detector.get_gene_biotype(left_gene, chrom=left_chr, pos=left_pos)
            right_biotype = self.detector.get_gene_biotype(right_gene, chrom=right_chr, pos=right_pos)
            if not self._is_allowed_biotype(left_biotype) or not self._is_allowed_biotype(right_biotype):
                fusions_to_discard.append(fusion_key)
        self._remove_discarded_fusions_internal(fusions_to_discard)

    def filter_multicopy_artifact_pairs(self) -> None:
        if not self.detector.fusion_metadata:
            return
        for _fusion_key, meta in self.detector.fusion_metadata.items():
            if meta.get("is_valid", True) is False:
                continue
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            both_multicopy = (
                left_gene
                and right_gene
                and self.is_multicopy_artifact_family(left_gene)
                and self.is_multicopy_artifact_family(right_gene)
            )
            both_ribo = (
                left_gene
                and right_gene
                and self._is_ribosomal_or_histone_gene(left_gene)
                and self._is_ribosomal_or_histone_gene(right_gene)
            )
            if both_multicopy or both_ribo:
                meta["is_valid"] = False
                meta.setdefault("reasons", []).append(
                    f"Both genes belong to multicopy artifact families: {left_gene} - {right_gene}"
                )

    def apply_frequency_filters(self) -> None:
        if not self.detector.fusion_metadata:
            return
        left_gene_count = defaultdict(int)
        right_gene_count = defaultdict(int)
        for meta in self.detector.fusion_metadata.values():
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            if left_gene and left_gene != "intergenic":
                left_gene_count[left_gene] += 1
            if right_gene and right_gene != "intergenic":
                right_gene_count[right_gene] += 1

        for meta in self.detector.fusion_metadata.values():
            if meta.get("is_valid", True) is False:
                continue
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            artifact_reasons = []
            for side, gene, counts in (("Left", left_gene, left_gene_count), ("Right", right_gene, right_gene_count)):
                if not gene or gene == "intergenic":
                    continue
                is_rib_hist = self._is_ribosomal_or_histone_gene(gene)
                threshold = 2 if is_rib_hist else 4
                count = counts.get(gene, 0)
                if count > threshold:
                    kind = "ribosomal/histone" if is_rib_hist else "other"
                    artifact_reasons.append(f"{side} gene '{gene}' appears {count} times (threshold={threshold} for {kind} genes)")
            if artifact_reasons:
                meta["is_valid"] = False
                meta.setdefault("reasons", []).extend(artifact_reasons)

    def _is_ribosomal_or_histone_gene(self, gene_name: Optional[str]) -> bool:
        if not gene_name or not isinstance(gene_name, str):
            return False
        return gene_name.upper().startswith(RIBOSOMAL_HISTONE_PREFIXES)

    def is_driver_gene(self, gene: Optional[str]) -> bool:
        return gene in DRIVER_GENES

    def is_multicopy_artifact_family(self, gene: str) -> bool:
        return gene.upper().startswith(MULTICOPY_PREFIXES)

    def confidence(self, meta: dict, flags: Optional[dict] = None) -> float:
        support = min(meta.get("support", 0), 10) / 10.0
        recon = 0.2 if meta.get("reconstruction_ok") else 0.0
        hits = meta.get("realignment_hits", 0)
        mapq = meta.get("best_hit_mapq", 0)
        realign = min(hits, 3) * 0.1 + min(mapq, 30) / 300.0
        priors = 0.0
        left = meta.get("left_gene")
        right = meta.get("right_gene")
        if left and right:
            if self.is_multicopy_artifact_family(left) or self.is_multicopy_artifact_family(right):
                priors -= 0.30
            if self.is_driver_gene(left) or self.is_driver_gene(right):
                priors += 0.20
        return max(0.0, min(1.0, support + recon + realign + priors))

    def _merge_fully_identical(self):
        bp_to_fusion_keys = defaultdict(list)
        for fusion_key, meta in self.detector.fusion_metadata.items():
            consensus_bp = meta.get("consensus_bp")
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            if not consensus_bp or not left_gene or not right_gene:
                continue
            bp_to_fusion_keys[(left_gene, right_gene, consensus_bp)].append(fusion_key)

        fusions_to_discard = set()
        for fusion_keys in bp_to_fusion_keys.values():
            if len(fusion_keys) < 2:
                continue
            keep_key = fusion_keys[0]
            for discard_key in fusion_keys[1:]:
                logger.info("Merging identical duplicates: %s <- %s", keep_key, discard_key)
                self._merge_fusion_candidates(keep_key, discard_key)
                fusions_to_discard.add(discard_key)
        return fusions_to_discard

    def _merge_fusion_candidates(self, keep_fusion_key: str, discard_fusion_key: str) -> None:
        if discard_fusion_key in self.detector.fusion_candidates:
            keep_reads = self.detector.fusion_candidates.get(keep_fusion_key, set())
            keep_reads.update(self.detector.fusion_candidates[discard_fusion_key])
            self.detector.fusion_candidates[keep_fusion_key] = keep_reads

        if discard_fusion_key in self.detector.fusion_breakpoints:
            keep_bps = self.detector.fusion_breakpoints.get(keep_fusion_key, {})
            for bp, count in self.detector.fusion_breakpoints[discard_fusion_key].items():
                keep_bps[bp] = keep_bps.get(bp, 0) + count
            self.detector.fusion_breakpoints[keep_fusion_key] = keep_bps

        if keep_fusion_key in self.detector.fusion_metadata and discard_fusion_key in self.detector.fusion_metadata:
            meta = self.detector.fusion_metadata[keep_fusion_key]
            discard_meta = self.detector.fusion_metadata[discard_fusion_key]
            if "supporting_reads" in meta and "supporting_reads" in discard_meta:
                meta["supporting_reads"].update(discard_meta["supporting_reads"])
                meta["support"] = len(meta["supporting_reads"])
