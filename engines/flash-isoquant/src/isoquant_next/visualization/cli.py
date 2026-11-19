"""Visualization CLI runner (`isoquant-next-visualize`) built on OutputConfig, Polars, and Agg."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Dict, List, Optional

from isoquant_next.core.error_codes import IsoQuantExitCode
from isoquant_next.core.logging_setup import configure_logging
from isoquant_next.visualization.dict_builder import DictionaryBuilder
from isoquant_next.visualization.output_config import OutputConfig
from isoquant_next.visualization.plotter import PlotOutput

logger = logging.getLogger("IsoQuant")


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="isoquant-next-visualize",
        description="Visualize IsoQuant / isoquant-next transcript models and isoform counts.",
    )
    parser.add_argument("--output", "-o", type=str, required=True, help="IsoQuant output directory")
    parser.add_argument("--visual_output", type=str, default=None, help="Directory for plots [<output>/visualization]")
    parser.add_argument("--gtf", type=str, default=None, help="Optional reference GTF/.iqidx override")
    parser.add_argument("--counts", action="store_true", default=False, help="Include count/abundance bar charts")
    parser.add_argument("--find_genes", action="store_true", default=False, help="Rank divergent genes across groups/samples")
    parser.add_argument("--num_genes", type=int, default=50, help="Maximum number of genes to plot [50]")
    parser.add_argument("--genes", nargs="+", default=None, help="Explicit gene IDs to visualize")
    parser.add_argument("--known_genes", type=str, default=None, help="File with gene IDs of interest (one per line)")
    parser.add_argument("--ref_only", action="store_true", default=False, help="Plot only reference transcripts")
    parser.add_argument("--format", type=str, default="pdf", choices=["pdf", "png", "svg"], help="Plot file format [pdf]")
    parser.add_argument("--debug", action="store_true", default=False)
    return parser.parse_args(argv)


def _load_known_genes(path: Optional[str]) -> Optional[List[str]]:
    if not path:
        return None
    p = Path(path)
    if not p.is_file():
        return None
    return [line.strip() for line in p.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]


def run_visualize(argv: List[str]) -> int:
    if not argv:
        parse_args(["--help"])
        return 0

    args = parse_args(argv)
    configure_logging(debug=args.debug)

    if not os.path.isdir(args.output):
        logger.critical("Output folder does not exist: %s", args.output)
        return int(IsoQuantExitCode.OUTPUT_NOT_FOUND)

    cfg = OutputConfig(
        args.output,
        use_counts=bool(args.counts or args.find_genes),
        ref_only=True if args.ref_only else None,
        gtf=args.gtf,
    )
    known = _load_known_genes(args.known_genes)
    target_genes: Optional[List[str]] = list(args.genes) if args.genes else known

    builder = DictionaryBuilder(cfg)
    gene_map = builder.build(gene_ids=target_genes)
    if not gene_map:
        logger.warning("No gene annotations found under %s", args.output)
        return int(IsoQuantExitCode.SUCCESS)

    tx_to_gene: Dict[str, str] = {
        tid: gid for gid, grec in gene_map.items() for tid in grec.transcripts
    }

    counts_by_tx: Dict[str, Dict[str, float]] = {}
    counts_source: Optional[Path] = None
    for s in cfg.samples:
        cpath = s.discovered_transcript_counts_tsv or s.transcript_counts_tsv
        if cpath is not None and cpath.is_file():
            counts_source = cpath
            if args.counts or args.find_genes:
                per_sample = DictionaryBuilder.load_transcript_counts(cpath)
                for tid, row in per_sample.items():
                    dest = counts_by_tx.setdefault(tid, {})
                    for k, v in row.items():
                        dest[f"{s.prefix}:{k}"] = float(v)

    selected_ids: List[str]
    if target_genes is not None:
        selected_ids = [g for g in target_genes if g in gene_map]
    elif args.find_genes and counts_source is not None:
        from isoquant_next.visualization.gene_ranking import rank_and_visualize_genes

        ranked = rank_and_visualize_genes(
            counts_source,
            tx_to_gene,
            num_genes=args.num_genes,
            known_genes=known,
        )
        selected_ids = [r.gene_id for r in ranked if r.gene_id in gene_map]
        if not selected_ids:
            selected_ids = sorted(gene_map.keys())[: args.num_genes]
    else:
        selected_ids = sorted(gene_map.keys())[: args.num_genes]

    out_plot_dir = Path(args.visual_output) if args.visual_output else Path(args.output) / "visualization"
    plotter = PlotOutput(out_plot_dir, image_format=args.format)
    plotter.plot_many(
        (gene_map[gid] for gid in selected_ids),
        counts_by_tx=counts_by_tx if (args.counts or args.find_genes) else None,
    )
    return int(IsoQuantExitCode.SUCCESS)
