"""Headless Agg Matplotlib gene model and isoform abundance plotter with leak-free `plt.close(fig)`."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=True)

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import matplotlib.patches as patches
import matplotlib.pyplot as plt

from isoquant_next.visualization.dict_builder import GeneModelRecord, TranscriptAnnotation


class PlotOutput:
    """Renders gene exon structures and isoform expression plots using non-interactive Agg backend."""

    def __init__(self, out_dir: str | Path, *, image_format: str = "pdf", dpi: int = 150) -> None:
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.image_format = image_format.lstrip(".").lower() or "pdf"
        self.dpi = int(dpi)

    @staticmethod
    def _sanitize_filename(name: str) -> str:
        return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name)

    def plot_gene_model(
        self,
        gene: GeneModelRecord,
        *,
        counts_by_tx: Optional[Dict[str, Dict[str, float]]] = None,
        sample_names: Optional[List[str]] = None,
    ) -> Optional[Path]:
        if not gene.transcripts:
            return None

        tx_list = sorted(
            gene.transcripts.values(),
            key=lambda t: (t.is_discovered, t.transcript_id),
        )
        n_tx = len(tx_list)
        g_start, g_end = gene.span
        if g_end <= g_start:
            g_end = g_start + 1

        fig_height = max(2.8, 1.2 + 0.45 * n_tx)
        has_counts = counts_by_tx is not None and len(counts_by_tx) > 0
        if has_counts:
            fig, (ax_model, ax_bar) = plt.subplots(
                1,
                2,
                figsize=(12, fig_height),
                gridspec_kw={"width_ratios": [3.2, 1.3]},
            )
        else:
            fig, ax_model = plt.subplots(1, 1, figsize=(10, fig_height))
            ax_bar = None

        try:
            pad = max(50.0, 0.04 * (g_end - g_start))
            ax_model.set_xlim(g_start - pad, g_end + pad)
            ax_model.set_ylim(-0.6, n_tx - 0.3)

            y_positions: List[float] = []
            y_labels: List[str] = []
            for row_idx, tx in enumerate(tx_list):
                y = float(n_tx - 1 - row_idx)
                y_positions.append(y)
                y_labels.append(tx.transcript_id)
                color = "#d95f02" if tx.is_discovered else "#1b9e77"
                if tx.exons:
                    ax_model.plot(
                        [tx.exons[0][0], tx.exons[-1][1]],
                        [y, y],
                        color=color,
                        linewidth=1.2,
                        zorder=1,
                    )
                    for ex_s, ex_e in tx.exons:
                        width = max(1, ex_e - ex_s + 1)
                        rect = patches.Rectangle(
                            (ex_s, y - 0.22),
                            width,
                            0.44,
                            facecolor=color,
                            edgecolor="#222222",
                            linewidth=0.5,
                            zorder=2,
                        )
                        ax_model.add_patch(rect)

            ax_model.set_yticks(y_positions)
            ax_model.set_yticklabels(y_labels, fontsize=8)
            ax_model.set_xlabel(f"{gene.chrom} ({gene.strand}) genomic coordinate")
            title_name = f"{gene.gene_name} ({gene.gene_id})" if gene.gene_name != gene.gene_id else gene.gene_id
            ax_model.set_title(title_name)
            ax_model.spines["top"].set_visible(False)
            ax_model.spines["right"].set_visible(False)

            if ax_bar is not None and counts_by_tx is not None:
                bar_vals: List[float] = []
                for tx in tx_list:
                    row = counts_by_tx.get(tx.transcript_id, {})
                    if sample_names:
                        val = sum(float(row.get(sn, 0.0)) for sn in sample_names)
                    else:
                        val = sum(float(v) for v in row.values())
                    bar_vals.append(val)
                bar_colors = ["#d95f02" if tx.is_discovered else "#1b9e77" for tx in tx_list]
                ax_bar.barh(y_positions, bar_vals, height=0.5, color=bar_colors, edgecolor="#333333", linewidth=0.4)
                ax_bar.set_ylim(-0.6, n_tx - 0.3)
                ax_bar.set_yticks(y_positions)
                ax_bar.set_yticklabels([])
                ax_bar.set_xlabel("Total counts / abundance")
                ax_bar.spines["top"].set_visible(False)
                ax_bar.spines["right"].set_visible(False)

            fig.tight_layout()
            out_path = self.out_dir / f"{self._sanitize_filename(gene.gene_id)}.{self.image_format}"
            fig.savefig(out_path, format=self.image_format, dpi=self.dpi)
            return out_path
        finally:
            plt.close(fig)

    def plot_count_summary_pie(self, stats_counts: Dict[str, float], *, filename: str = "summary_pie") -> Optional[Path]:
        """Draw non-interactive assignment/transcript category summary chart and close figure cleanly."""
        filtered = [(k, float(v)) for k, v in stats_counts.items() if float(v) > 0]
        if not filtered:
            return None
        labels = [k for k, _ in filtered]
        values = [v for _, v in filtered]
        fig, ax = plt.subplots(figsize=(6, 5))
        try:
            ax.pie(values, labels=labels, autopct="%1.1f%%", startangle=90)
            ax.set_title("Feature / Assignment Breakdown")
            out_path = self.out_dir / f"{self._sanitize_filename(filename)}.{self.image_format}"
            fig.savefig(out_path, format=self.image_format, dpi=self.dpi)
            return out_path
        finally:
            plt.close(fig)

    def plot_many(
        self,
        genes: Iterable[GeneModelRecord],
        *,
        counts_by_tx: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> List[Path]:
        written: List[Path] = []
        for g in genes:
            p = self.plot_gene_model(g, counts_by_tx=counts_by_tx)
            if p is not None:
                written.append(p)
        return written
