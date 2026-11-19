"""Stage 5 of transcript discovery: groups overlapping novel/known models into genes."""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, Tuple

from isoquant_next.annotation.gene_info import TranscriptModelType
from isoquant_next.core.utils import intersection_len, interval_len, junctions_from_blocks

logger = logging.getLogger("IsoQuant")


class TranscriptToGeneJoiner:
    """Merges same-strand genes sharing splice junctions or strongly overlapping coordinates."""

    STRONG_OVERLAP_FRACTION = 0.5

    def __init__(self, transcript_model_storage, gene_info) -> None:
        self.gene_info = gene_info
        self.transcript_model_storage = transcript_model_storage
        self.gene_introns: Dict[str, set] = defaultdict(set)
        self.gene_strands: Dict[str, str] = {}
        self.gene_regions: Dict[str, Tuple[int, int]] = {}
        self.gene_to_transcripts: Dict[str, set[str]] = defaultdict(set)

        ref_regions = self.gene_info.get_gene_regions()
        for gene_id, strand in self.gene_info.gene_strands.items():
            self.gene_strands[gene_id] = strand
            self.gene_regions[gene_id] = ref_regions[gene_id]
        for tid, gid in self.gene_info.gene_id_map.items():
            self.gene_introns[gid].update(self.gene_info.all_isoforms_introns[tid])
            self.gene_to_transcripts[gid].add(tid)

        for t in self.transcript_model_storage:
            if t.transcript_type == TranscriptModelType.known:
                continue
            region = (t.exon_blocks[0][0], t.exon_blocks[-1][1])
            if t.gene_id not in self.gene_regions:
                self.gene_regions[t.gene_id] = region
                self.gene_strands[t.gene_id] = t.strand
            else:
                prev = self.gene_regions[t.gene_id]
                self.gene_regions[t.gene_id] = (min(prev[0], region[0]), max(prev[1], region[1]))
                assert self.gene_strands[t.gene_id] == t.strand
            self.gene_introns[t.gene_id].update(junctions_from_blocks(t.exon_blocks))
            self.gene_to_transcripts[t.gene_id].add(t.transcript_id)
        self.scores: Dict[Tuple[str, str], float] = {}

    def count_score(self, gene1: str, gene2: str) -> float:
        if self.gene_strands[gene1] != self.gene_strands[gene2]:
            return 0.0
        introns1 = self.gene_introns[gene1]
        introns2 = self.gene_introns[gene2]
        shared = len(introns1 & introns2)
        if shared > 0:
            return 1.0 + shared / max(1, len(introns1 | introns2))
        r1 = self.gene_regions[gene1]
        r2 = self.gene_regions[gene2]
        inter = intersection_len(r1, r2)
        if inter <= 0:
            return 0.0
        overlap_frac = inter / min(interval_len(r1), interval_len(r2))
        return overlap_frac if overlap_frac >= self.STRONG_OVERLAP_FRACTION else 0.0

    def count_scores(self) -> None:
        gids = list(self.gene_to_transcripts.keys())
        for g1 in gids:
            for g2 in gids:
                if g1 == g2 or (g1 in self.gene_info.gene_strands and g2 in self.gene_info.gene_strands):
                    continue
                pair = tuple(sorted((g1, g2)))
                if pair not in self.scores:
                    self.scores[pair] = self.count_score(g1, g2)

    def merge_genes(self, gene1: str, gene2: str) -> None:
        self.gene_regions[gene1] = (
            min(self.gene_regions[gene1][0], self.gene_regions[gene2][0]),
            max(self.gene_regions[gene1][1], self.gene_regions[gene2][1]),
        )
        self.gene_introns[gene1].update(self.gene_introns[gene2])
        self.gene_to_transcripts[gene1].update(self.gene_to_transcripts[gene2])
        del self.gene_regions[gene2]
        del self.gene_introns[gene2]
        del self.gene_to_transcripts[gene2]
        del self.gene_strands[gene2]

        new_scores: Dict[Tuple[str, str], float] = {}
        for pair, val in self.scores.items():
            if gene2 in pair:
                continue
            new_scores[pair] = self.count_score(*pair) if gene1 in pair else val
        self.scores = new_scores

    def join_transcripts(self):
        self.count_scores()
        while len(self.scores) > 1:
            best_pair = max(self.scores, key=self.scores.get)
            if self.scores[best_pair] <= 0.0:
                break
            if best_pair[0] in self.gene_info.gene_strands:
                assert best_pair[1] not in self.gene_info.gene_strands
                self.merge_genes(best_pair[0], best_pair[1])
            else:
                assert best_pair[0] not in self.gene_info.gene_strands
                self.merge_genes(best_pair[1], best_pair[0])

        tx_to_gene: Dict[str, str] = {}
        for gene_id, tids in self.gene_to_transcripts.items():
            for tid in tids:
                tx_to_gene[tid] = gene_id
        for model in self.transcript_model_storage:
            model.gene_id = tx_to_gene[model.transcript_id]
        return self.transcript_model_storage
