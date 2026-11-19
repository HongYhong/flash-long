"""IsoQuant-next visualization package (`run_manifest.json`, `.iqidx`, Polars, Agg Matplotlib)."""

from __future__ import annotations

from typing import Any

from isoquant_next.visualization.cli import run_visualize
from isoquant_next.visualization.dict_builder import DictionaryBuilder, GeneModelRecord, TranscriptAnnotation
from isoquant_next.visualization.output_config import OutputConfig, SamplePaths
from isoquant_next.visualization.plotter import PlotOutput

__all__ = [
    "DictionaryBuilder",
    "GeneModelRecord",
    "OutputConfig",
    "PlotOutput",
    "RankedGene",
    "SamplePaths",
    "TranscriptAnnotation",
    "rank_and_visualize_genes",
    "run_visualize",
    "simplify_and_sum_transcripts",
]


def __getattr__(name: str) -> Any:
    if name in ("RankedGene", "rank_and_visualize_genes", "simplify_and_sum_transcripts"):
        from isoquant_next.visualization import gene_ranking

        return getattr(gene_ranking, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
