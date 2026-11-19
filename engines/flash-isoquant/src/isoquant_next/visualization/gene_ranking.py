"""Polars / NumPy differential isoform gene ranking without Pandas CoW pitfalls."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import polars as pl
from scipy.spatial.distance import euclidean


@dataclass(frozen=True)
class RankedGene:
    gene_id: str
    divergence_score: float
    total_counts: float
    num_isoforms: int


def _load_count_dataframe(path: str | Path) -> pl.DataFrame:
    df = pl.read_csv(
        path,
        separator="\t",
        comment_prefix=None,
        infer_schema_length=256,
    )
    first_col = df.columns[0]
    clean_first = first_col.lstrip("#")
    if clean_first != first_col:
        df = df.rename({first_col: clean_first})
    return df.filter(~pl.col(clean_first).cast(pl.Utf8).str.starts_with("__"))


def simplify_and_sum_transcripts(
    count_table: pl.DataFrame,
    transcript_to_gene: Dict[str, str],
    *,
    min_total_count: float = 1.0,
) -> Tuple[pl.DataFrame, Dict[str, np.ndarray]]:
    """Aggregate normalized isoform expression vectors per gene using Polars (zero Pandas)."""
    id_col = count_table.columns[0]
    sample_cols = [c for c in count_table.columns[1:]]
    if not sample_cols:
        return pl.DataFrame({id_col: []}), {}

    # Map transcript_id -> gene_id and cast sample counts cleanly as float64
    tx_ids = count_table[id_col].cast(pl.Utf8).to_list()
    gene_ids = [transcript_to_gene.get(tid, "") for tid in tx_ids]

    cast_exprs = [pl.col(c).cast(pl.Float64, strict=False).fill_null(0.0) for c in sample_cols]
    enriched = (
        count_table.with_columns(pl.Series("gene_id", gene_ids), *cast_exprs)
        .filter(pl.col("gene_id") != "")
    )
    if enriched.height == 0:
        return enriched.select(["gene_id", *sample_cols]), {}

    gene_totals = enriched.group_by("gene_id").agg([pl.col(c).sum() for c in sample_cols])
    # Filter genes meeting minimum read threshold
    totals_arr = gene_totals.select(sample_cols).to_numpy()
    keep_mask = np.asarray(totals_arr.sum(axis=1) >= min_total_count, dtype=bool)
    kept_genes = set(g for g, keep in zip(gene_totals["gene_id"].to_list(), keep_mask) if keep)

    filtered = enriched.filter(pl.col("gene_id").is_in(list(kept_genes)))
    per_gene_mats: Dict[str, np.ndarray] = {}
    for gname, sub in filtered.group_by(["gene_id"], maintain_order=True):
        gid = gname[0] if isinstance(gname, tuple) else gname
        per_gene_mats[str(gid)] = np.asarray(sub.select(sample_cols).to_numpy(), dtype=np.float64)
    return gene_totals, per_gene_mats


def rank_and_visualize_genes(
    counts_path: str | Path,
    transcript_to_gene: Dict[str, str],
    *,
    num_genes: int = 100,
    known_genes: Optional[Iterable[str]] = None,
    min_total_count: float = 5.0,
) -> List[RankedGene]:
    """Rank genes by cross-condition normalized isoform proportion divergence (Euclidean)."""
    df = _load_count_dataframe(counts_path)
    _gene_totals, per_gene_mats = simplify_and_sum_transcripts(
        df, transcript_to_gene, min_total_count=min_total_count
    )
    allowed = set(known_genes) if known_genes is not None else None

    scored: List[RankedGene] = []
    for gid, mat in per_gene_mats.items():
        if allowed is not None and gid not in allowed:
            continue
        if mat.shape[0] < 2 or mat.shape[1] < 2:
            total = float(mat.sum())
            scored.append(RankedGene(gene_id=gid, divergence_score=0.0, total_counts=total, num_isoforms=mat.shape[0]))
            continue
        col_sums = np.maximum(mat.sum(axis=0, keepdims=True), 1e-12)
        props = mat / col_sums
        # Pairwise max Euclidean divergence across columns
        max_div = 0.0
        ncols = props.shape[1]
        for i in range(ncols):
            for j in range(i + 1, ncols):
                dist = float(euclidean(props[:, i], props[:, j]))
                if dist > max_div:
                    max_div = dist
        scored.append(
            RankedGene(
                gene_id=gid,
                divergence_score=max_div,
                total_counts=float(mat.sum()),
                num_isoforms=int(mat.shape[0]),
            )
        )

    scored.sort(key=lambda r: (-r.divergence_score, -r.total_counts, r.gene_id))
    return scored[: max(0, int(num_genes))]
