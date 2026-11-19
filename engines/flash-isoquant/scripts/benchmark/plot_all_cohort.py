#!/usr/bin/env python3
"""Generate publication-quality full-cohort benchmark, memory, correlation, and feature figures."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np


C_ORIG = "#64748B"
C_NEXT = "#0284C7"
C_ACCENT = "#059669"
C_WARN = "#D97706"
C_PURPLE = "#7C3AED"
C_INDIGO = "#1D4ED8"

PROTOCOL_ORDER = [
    ("cDNA_PacBio", "cDNA PacBio"),
    ("cDNA_PacBio_v46", "cDNA PacBio (v46)"),
    ("cDNA_PacBio_vM36", "cDNA PacBio (vM36)"),
    ("capped_PacBio", "Capped PacBio"),
    ("directRNA", "directRNA (ONT)"),
    ("cDNA_ONT", "cDNA ONT"),
    ("capped_ONT", "Capped ONT"),
    ("R2C2_ONT", "R2C2 ONT"),
    ("R2C2_ONT_rmdup", "R2C2 ONT (rmdup)"),
    ("RBMX_KD", "RBMX KD (v46)"),
    ("SIRV_simulation", "SIRV Simulation (v46)"),
    ("Sequin_SIRV4", "Sequin + SIRV4"),
    ("MultiSample_YAML", "Multi-Sample YAML"),
    ("FASTQ_directRNA", "FASTQ -> minimap2"),
    ("Illumina_Hybrid", "Illumina Hybrid"),
]


def _fmt_duration(sec: float) -> str:
    if sec >= 3600:
        return f"{sec / 3600:.2f}h"
    if sec >= 120:
        return f"{sec / 60:.1f}m"
    return f"{sec:.1f}s"


def _load_counts(path: Path) -> Dict[str, float]:
    out: Dict[str, float] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#") or line.startswith("__"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                try:
                    out[p[0]] = float(p[1])
                except ValueError:
                    pass
    return out


def _load_pss_trace(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    t, rss, pss = [], [], []
    with open(path, "r", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                t.append(float(p[0]) / 60.0)
                rss.append(int(p[2]) / 1048576.0)
                pss.append(int(p[3]) / 1048576.0)
    return np.array(t), np.array(rss), np.array(pss)


def _group_paired_by_protocol(rows: List[dict]) -> List[Tuple[str, str, List[dict]]]:
    by_p: Dict[str, List[dict]] = defaultdict(list)
    for r in rows:
        if r.get("is_paired_engine_cmp"):
            by_p[r["protocol"]].append(r)
    out = []
    for pkey, plabel in PROTOCOL_ORDER:
        if by_p.get(pkey):
            out.append((pkey, plabel, by_p[pkey]))
    return out


def plot_fig1_speedup(cohort: dict, out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.3, "axes.edgecolor": "#CBD5E1"})
    groups = _group_paired_by_protocol(cohort["rows"])
    n = len(groups)
    y = np.arange(n)[::-1]
    h = 0.34
    rng = np.random.default_rng(42)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.0, max(7.8, 0.52 * n + 2.2)), gridspec_kw={"width_ratios": [1.3, 1.0]}, sharey=True)

    orig_meds, next_meds, spd_meds, labels = [], [], [], []
    all_spds: List[float] = []
    for yi, (_, plabel, items) in zip(y, groups):
        ow = np.array([x["orig_wall_s"] for x in items])
        nw = np.array([x["next_wall_s"] for x in items])
        sp = np.array([x["speedup_x"] for x in items])
        all_spds.extend(sp.tolist())
        orig_meds.append(float(np.median(ow)))
        next_meds.append(float(np.median(nw)))
        spd_meds.append(float(np.median(sp)))
        labels.append(f"{plabel} (n={len(items)})")

        # Jittered replicate points
        jit_o = rng.uniform(-0.10, 0.10, size=len(ow))
        jit_n = rng.uniform(-0.10, 0.10, size=len(nw))
        ax1.scatter(ow, yi + h / 2 + jit_o, s=16, color="#1E293B", alpha=0.55, zorder=4)
        ax1.scatter(nw, yi - h / 2 + jit_n, s=16, color="#0369A1", alpha=0.75, zorder=4)
        ax2.scatter(sp, yi + jit_n, s=20, color="#065F46" if np.median(sp) >= 10 else "#0369A1", alpha=0.72, zorder=4)

        # Min-max range whiskers on both runtime and speedup
        if len(ow) > 1:
            ax1.plot([ow.min(), ow.max()], [yi + h / 2, yi + h / 2], color="#334155", linewidth=1.1, alpha=0.55, zorder=3)
            ax1.plot([nw.min(), nw.max()], [yi - h / 2, yi - h / 2], color="#0369A1", linewidth=1.1, alpha=0.65, zorder=3)
        if len(sp) > 1:
            ax2.plot([sp.min(), sp.max()], [yi, yi], color="#0F172A", linewidth=1.3, alpha=0.6, zorder=3)

    ax1.barh(y + h / 2, orig_meds, height=h, color=C_ORIG, alpha=0.88, label="Original isoquant.py (Median + Replicates)")
    ax1.barh(y - h / 2, next_meds, height=h, color=C_NEXT, alpha=0.92, label="Refactored isoquant-next (Median + Replicates)")
    ax1.set_xscale("log")
    max_ow = max(r["orig_wall_s"] for _, _, items in groups for r in items)
    min_nw = min(r["next_wall_s"] for _, _, items in groups for r in items)
    ax1.set_xlim(max(1.0, min_nw * 0.55), max_ow * 4.5)
    ax1.set_ylim(-0.7, n - 0.3)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontweight="medium")
    ax1.set_xlabel("Wall-Clock Runtime (seconds, log scale)")
    ax1.set_title("A. End-to-End Runtime Across Full Slurm Cohort", fontweight="bold", pad=10)
    ax1.grid(axis="x", linestyle="--", alpha=0.38)

    for yi, (_, _, items), om, nm in zip(y, groups, orig_meds, next_meds):
        omax = max(x["orig_wall_s"] for x in items)
        nmax = max(x["next_wall_s"] for x in items)
        ax1.text(max(om, omax) * 1.15, yi + h / 2, _fmt_duration(om), va="center", ha="left", fontsize=7.5, color="#334155")
        ax1.text(max(nm, nmax) * 1.15, yi - h / 2, _fmt_duration(nm), va="center", ha="left", fontsize=7.5, color="#0369A1", fontweight="bold")

    bar_cols = [C_ACCENT if s >= 10.0 else C_NEXT for s in spd_meds]
    ax2.barh(y, spd_meds, height=0.54, color=bar_cols, alpha=0.85)
    ax2.axvline(1.0, color="#94A3B8", linestyle="--", linewidth=1.1)
    max_sp = max(all_spds)
    ax2.set_xlim(0, max_sp * 1.24)
    ax2.set_xlabel("Speedup Factor (x faster than isoquant.py)")
    ax2.set_title(f"B. Measured Speedup Distribution ({len(all_spds)} Paired Runs)", fontweight="bold", pad=10)
    ax2.grid(axis="x", linestyle="--", alpha=0.38)
    for yi, (_, _, items), sm in zip(y, groups, spd_meds):
        sp = [x["speedup_x"] for x in items]
        hi = max(sp)
        txt = f"{sm:.2f}x" if len(sp) == 1 else f"{sm:.2f}x [{min(sp):.1f}–{hi:.1f}x]"
        ax2.text(max(sm, hi) + max_sp * 0.02, yi, txt, va="center", ha="left", fontsize=7.8, fontweight="bold", color="#0F172A")

    handles, lbls = ax1.get_legend_handles_labels()
    fig.legend(handles, lbls, loc="upper center", bbox_to_anchor=(0.56, 0.942), ncol=2, frameon=True, framealpha=0.95)
    g = cohort["global"]
    fig.suptitle(
        f"IsoQuant Scheme B Refactoring: Full-Cohort Runtime & Speedup ({g['paired_engine_tasks']} Paired Slurm Runs, Median {g['overall_speedup_median_x']}x)",
        fontsize=12.8,
        fontweight="bold",
        y=0.985,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_fig2_memory_and_scaling(cohort: dict, out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.1, "axes.edgecolor": "#CBD5E1"})
    groups = _group_paired_by_protocol(cohort["rows"])
    fig = plt.figure(figsize=(16.0, 9.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.35, 1.0])
    ax_pss = fig.add_subplot(gs[0, :])
    ax_rss = fig.add_subplot(gs[1, 0])
    ax_scale = fig.add_subplot(gs[1, 1])

    x = np.arange(len(groups))
    w = 0.36
    opss = [float(np.median([r["orig_pss_gib"] for r in items])) for _, _, items in groups]
    npss = [float(np.median([r["next_pss_gib"] for r in items])) for _, _, items in groups]
    short_names = [f"{lbl.replace(' -> ', chr(10)).replace(' ', chr(10))}\n(n={len(items)})" for _, lbl, items in groups]

    ax_pss.bar(x - w / 2, opss, width=w, color=C_ORIG, label="Original isoquant.py Median Peak Tree PSS (GiB)")
    ax_pss.bar(x + w / 2, npss, width=w, color=C_NEXT, label="Refactored isoquant-next Median Peak Tree PSS (GiB)")
    rng = np.random.default_rng(7)
    for xi, (_, _, items) in zip(x, groups):
        ov = [r["orig_pss_gib"] for r in items]
        nv = [r["next_pss_gib"] for r in items]
        ax_pss.scatter(np.full(len(ov), xi - w / 2) + rng.uniform(-0.07, 0.07, len(ov)), ov, s=11, color="#1E293B", alpha=0.5, zorder=4)
        ax_pss.scatter(np.full(len(nv), xi + w / 2) + rng.uniform(-0.07, 0.07, len(nv)), nv, s=11, color="#0369A1", alpha=0.7, zorder=4)

    ax_pss.set_xticks(x)
    ax_pss.set_xticklabels(short_names, fontsize=7.3)
    ax_pss.set_ylabel("Peak Process-Tree PSS (GiB)")
    ax_pss.set_ylim(0, max(max(opss), max(npss)) * 1.25)
    ax_pss.set_title("A. Peak Process-Tree Memory (PSS via /proc/<pid>/smaps_rollup) Across All Protocols", fontweight="bold")
    ax_pss.grid(axis="y", linestyle="--", alpha=0.38)
    ax_pss.legend(loc="upper left", frameon=True)

    for xi, op, np_ in zip(x, opss, npss):
        pct = (op - np_) / op * 100.0 if op > 0 else 0.0
        txt = f"-{pct:.0f}%" if pct >= 1.0 else (f"+{-pct:.1f}%" if pct < -1.0 else "~0%")
        col = C_ACCENT if pct >= 50.0 else ("#0369A1" if pct > 5.0 else "#64748B")
        ax_pss.text(
            xi + w / 2, np_ + 0.65, txt, ha="center", va="bottom", fontsize=7.7, fontweight="bold", color=col,
            bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="none", alpha=0.85), zorder=5,
        )

    # Panel B: High-concurrency RSS comparison across major heavy cohorts
    rss_protos = [
        ("cDNA_PacBio", "cDNA\nPacBio"),
        ("directRNA", "directRNA\n(ONT)"),
        ("cDNA_ONT", "cDNA\nONT"),
        ("capped_ONT", "Capped\nONT"),
        ("R2C2_ONT", "R2C2\nONT"),
        ("RBMX_KD", "RBMX KD\n(v46)"),
        ("SIRV_simulation", "SIRV Sim\n(v46)"),
    ]
    by_p = {k: items for k, _, items in groups}
    sel_rss = [(lbl, by_p[k]) for k, lbl in rss_protos if k in by_p]
    xr = np.arange(len(sel_rss))
    orss = [float(np.median([r["orig_rss_gib"] for r in items])) for _, items in sel_rss]
    nrss = [float(np.median([r["next_rss_gib"] for r in items])) for _, items in sel_rss]
    ax_rss.bar(xr - w / 2, orss, width=w, color=C_ORIG, label="Original Median Peak Tree RSS (GiB)")
    ax_rss.bar(xr + w / 2, nrss, width=w, color=C_NEXT, label="Refactored Median Peak Tree RSS (GiB)")
    ax_rss.set_xticks(xr)
    ax_rss.set_xticklabels([lbl for lbl, _ in sel_rss], fontsize=8.2)
    ax_rss.set_ylabel("Median Peak Tree RSS (GiB)")
    ax_rss.set_ylim(0, max(orss) * 1.22)
    ax_rss.set_title("B. High-Concurrency RSS Savings (Shared mmap .iqidx vs Fork Copy-on-Write)", fontweight="bold")
    ax_rss.grid(axis="y", linestyle="--", alpha=0.38)
    ax_rss.legend(loc="upper left", frameon=True, fontsize=8.0)
    for xi, orv, nrv in zip(xr, orss, nrss):
        red = (orv - nrv) / orv * 100.0
        ax_rss.text(xi - w / 2, orv + 0.6, f"{orv:.1f}G", ha="center", va="bottom", fontsize=7.6, color="#334155")
        ax_rss.text(xi + w / 2, nrv + 0.6, f"{nrv:.1f}G\n(-{red:.0f}%)", ha="center", va="bottom", fontsize=7.6, fontweight="bold", color=C_ACCENT)

    # Panel C: ENCFF400BQQ thread scaling + warm .iqidx
    by_slug = {r["slug"]: r for r in cohort["rows"]}
    req_scale = ["g1_encff400bqq_t1", "g1_encff400bqq_t8_abba3", "g1_encff400bqq_t16"]
    if all(s in by_slug for s in req_scale):
        threads = [1, 8, 16]
        o_times = [by_slug[s]["orig_wall_s"] / 60.0 for s in req_scale]
        n_times = [by_slug[s]["next_wall_s"] / 60.0 for s in req_scale]
        ax_scale.plot(threads, o_times, "o-", color=C_ORIG, linewidth=2.2, markersize=7, label="Original isoquant.py (Cold GTF -> DB)")
        ax_scale.plot(threads, n_times, "s-", color=C_NEXT, linewidth=2.2, markersize=7, label="Refactored isoquant-next (Cold GTF -> .iqidx)")
        if "g1_encff400bqq_t8_warm_abba3" in by_slug:
            warm_min = by_slug["g1_encff400bqq_t8_warm_abba3"]["next_wall_s"] / 60.0
            warm_spd = by_slug["g1_encff400bqq_t8_warm_abba3"]["speedup_x"]
            ax_scale.scatter([8], [warm_min], marker="*", s=200, color=C_ACCENT, zorder=5, label=f"Warm .iqidx (-t 8: {warm_min:.2f}m, {warm_spd:.1f}x)")
            ax_scale.text(8, warm_min * 0.72, f"Warm .iqidx: {warm_min:.2f}m ({warm_spd:.1f}x)", ha="center", va="top", fontsize=8.0, fontweight="bold", color=C_ACCENT)
        ax_scale.set_yscale("log")
        ax_scale.set_ylim(0.55, 280)
        ax_scale.set_xlim(-1.0, 17.8)
        ax_scale.set_xticks(threads)
        ax_scale.set_xticklabels(["-t 1", "-t 8", "-t 16"])
        ax_scale.set_ylabel("Wall-Clock Runtime (minutes, log scale)")
        ax_scale.set_title("C. Thread Scaling on ENCFF400BQQ (-t 1 -> 8 -> 16 + Warm .iqidx)", fontweight="bold")
        ax_scale.grid(True, which="both", linestyle="--", alpha=0.35)
        ax_scale.legend(loc="upper right", frameon=True, fontsize=7.9)
        for th, ot, nt in zip(threads, o_times, n_times):
            ha = "left" if th == 1 else "center"
            xt = th - 0.1 if th == 1 else th
            ax_scale.text(xt, ot * 1.18, f"{ot:.1f}m", ha=ha, va="bottom", fontsize=8.1, color="#334155")
            ax_scale.text(xt, nt * 1.24, f"{nt:.2f}m ({ot/nt:.1f}x)", ha=ha, va="bottom", fontsize=7.9, fontweight="bold", color="#0369A1")

    fig.suptitle("Full-Cohort Process-Tree Memory (PSS / RSS) & Thread Scaling", fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_fig3_correlation_and_sirv(cohort: dict, scatter_dir: Path | None, bench_root: Path, out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.1, "axes.edgecolor": "#CBD5E1"})
    fig, axes = plt.subplots(2, 2, figsize=(15.8, 11.4))
    ax_g, ax_t, ax_par, ax_sim = axes.ravel()

    # Panels A & B: Gene and Transcript raw count parity scatter if scatter_dir provided, else summary correlation vs speedup
    has_scatter = scatter_dir is not None and (scatter_dir / "orig_gene.tsv").exists() and (scatter_dir / "next_gene.tsv").exists()
    if has_scatter:
        for ax, fname_o, fname_n, kind, color in (
            (ax_g, "orig_gene.tsv", "next_gene.tsv", "Gene", C_NEXT),
            (ax_t, "orig_tx.tsv", "next_tx.tsv", "Transcript", C_ACCENT),
        ):
            ov = _load_counts(scatter_dir / fname_o)
            nv = _load_counts(scatter_dir / fname_n)
            keys = sorted(set(ov) | set(nv))
            xs = np.array([ov.get(k, 0.0) for k in keys])
            ys = np.array([nv.get(k, 0.0) for k in keys])
            mask = (xs > 0) | (ys > 0)
            r = float(np.corrcoef(xs, ys)[0, 1]) if len(xs) > 1 else 1.0
            exact = float(np.mean(np.abs(xs - ys) <= 1e-6))
            lx, ly = np.log10(xs[mask] + 1.0), np.log10(ys[mask] + 1.0)
            lim = max(lx.max(), ly.max()) * 1.06 if len(lx) else 1.0
            ax.scatter(lx, ly, s=11, alpha=0.45, color=color, edgecolors="none", label=f"Expressed {kind.lower()}s (n={mask.sum():,})")
            ax.plot([0, lim], [0, lim], "--", color="#0F172A", linewidth=1.1, label="Identity y = x")
            ax.set_xlim(-0.08, lim)
            ax.set_ylim(-0.08, lim)
            ax.set_xlabel(f"Original isoquant.py: log10({kind} Count + 1)")
            ax.set_ylabel(f"Refactored isoquant-next: log10({kind} Count + 1)")
            ax.set_title(f"{'A' if kind=='Gene' else 'B'}. Full-Genome {kind} Quantification Correlation", fontweight="bold")
            ax.grid(True, linestyle="--", alpha=0.35)
            ax.text(
                0.04, 0.88,
                f"Pearson r (raw) = {r:.6f}\nExact match (all {len(keys):,}) = {exact*100:.2f}%\nMax |diff| = {np.max(np.abs(xs - ys)):.4g}",
                transform=ax.transAxes, va="top",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8FAFC", edgecolor="#CBD5E1"),
                fontsize=8.4,
            )
            ax.legend(loc="lower right", frameon=True)
    else:
        paired = [r for r in cohort["rows"] if r.get("is_paired_engine_cmp")]
        g_vals = [r["primary_gene_or_disc_r"] for r in paired if r.get("primary_gene_or_disc_r") is not None]
        t_vals = [r["primary_tx_or_disc_r"] for r in paired if r.get("primary_tx_or_disc_r") is not None]
        ax_g.hist(g_vals, bins=15, color=C_NEXT, alpha=0.85, edgecolor="white")
        ax_g.set_title("A. Cohort Gene / Discovered Pearson r Distribution", fontweight="bold")
        ax_g.set_xlabel("Pearson r")
        ax_t.hist(t_vals, bins=15, color=C_ACCENT, alpha=0.85, edgecolor="white")
        ax_t.set_title("B. Cohort Transcript Pearson r Distribution", fontweight="bold")
        ax_t.set_xlabel("Pearson r")

    # Panel C: Quantitative Parity Metrics Grouped by Protocol
    groups = _group_paired_by_protocol(cohort["rows"])
    yp = np.arange(len(groups))[::-1]
    rng = np.random.default_rng(11)
    y_labels = []
    min_tr, min_ri, min_ea, min_tm = 1.0, 1.0, 1.0, 1.0
    for yi, (_, plabel, items) in zip(yp, groups):
        y_labels.append(f"{plabel} (n={len(items)})")
        tr_vals = [r["primary_tx_or_disc_r"] for r in items if r.get("primary_tx_or_disc_r") is not None]
        ri_vals = [r["read_info_full_agreement"] for r in items if r.get("read_info_full_agreement") is not None]
        ea_vals = [r["extended_annotation_jaccard"] for r in items if r.get("extended_annotation_jaccard") is not None]
        tm_vals = [r["transcript_models_jaccard"] for r in items if r.get("transcript_models_jaccard") is not None]
        if tr_vals:
            min_tr = min(min_tr, min(tr_vals))
            ax_par.scatter(tr_vals, np.full(len(tr_vals), yi) + rng.uniform(-0.12, 0.12, len(tr_vals)), s=24, color=C_INDIGO, marker="o", alpha=0.8)
        if ri_vals:
            min_ri = min(min_ri, min(ri_vals))
            ax_par.scatter(ri_vals, np.full(len(ri_vals), yi) + rng.uniform(-0.12, 0.12, len(ri_vals)), s=22, color=C_ACCENT, marker="s", alpha=0.75)
        if ea_vals:
            min_ea = min(min_ea, min(ea_vals))
            ax_par.scatter(ea_vals, np.full(len(ea_vals), yi) + rng.uniform(-0.12, 0.12, len(ea_vals)), s=20, color=C_PURPLE, marker="D", alpha=0.75)
        if tm_vals:
            min_tm = min(min_tm, min(tm_vals))
            ax_par.scatter(tm_vals, np.full(len(tm_vals), yi) + rng.uniform(-0.12, 0.12, len(tm_vals)), s=22, color=C_WARN, marker="^", alpha=0.75)

    # Dummy handles for clean legend
    ax_par.scatter([], [], marker="o", s=38, color=C_INDIGO, label=f"Tx / Disc. Pearson r (min {min_tr:.6f})")
    ax_par.scatter([], [], marker="s", s=36, color=C_ACCENT, label=f"read_info Agreement (min {min_ri:.6f})")
    ax_par.scatter([], [], marker="D", s=34, color=C_PURPLE, label=f"Extended GTF Jaccard (min {min_ea:.6f})")
    ax_par.scatter([], [], marker="^", s=38, color=C_WARN, label=f"Transcript Models Jaccard (min {min_tm:.6f})")
    ax_par.axvline(1.0, color="#0F172A", linestyle=":", linewidth=1.0)
    low_x = max(0.90, min(min_tr, min_ri, min_ea, min_tm) - 0.004)
    ax_par.set_xlim(low_x, 1.0005)
    ax_par.set_yticks(yp)
    ax_par.set_yticklabels(y_labels, fontsize=7.6)
    ax_par.set_xlabel(f"Correlation / Agreement / Structural Jaccard ({low_x:.4f} – 1.0000)")
    ax_par.set_title("C. Cohort-Wide Output Parity Across All Protocols", fontweight="bold")
    ax_par.grid(True, linestyle="--", alpha=0.35)
    ax_par.legend(loc="lower left", frameon=True, framealpha=0.95, fontsize=7.4)

    # Panel D: All available SIRV_simulation replicates (simulated1..3)
    sim_docs = []
    for i in (1, 2, 3):
        ep = bench_root / f"g3_simulated{i}_t16" / f"eval_simulated{i}.json"
        if ep.is_file():
            sim_docs.append(json.loads(ep.read_text(encoding="utf-8")))
    metrics_keys = [
        ("read_precision", "Read\nPrecision"),
        ("read_recall", "Read\nRecall"),
        ("model_known_over_all_precision", "Model\nPrecision"),
        ("model_recall_ge2_reads", "Model\nRecall(>=2r)"),
        ("quant_precision", "Quant\nPrecision"),
        ("quant_recall_ge2_reads", "Quant\nRecall(>=2r)"),
    ]
    if sim_docs:
        xs = np.arange(len(metrics_keys))
        bw = 0.36
        b_vals = [float(np.mean([d["baseline"][k] for d in sim_docs])) for k, _ in metrics_keys]
        c_vals = [float(np.mean([d["candidate"][k] for d in sim_docs])) for k, _ in metrics_keys]
        ax_sim.bar(xs - bw / 2, b_vals, width=bw, color=C_ORIG, label=f"Original isoquant.py (mean n={len(sim_docs)})")
        ax_sim.bar(xs + bw / 2, c_vals, width=bw, color=C_NEXT, label=f"Refactored isoquant-next (mean n={len(sim_docs)})")
        for rep_i, d in enumerate(sim_docs):
            ax_sim.scatter(xs - bw / 2 + (rep_i - 1) * 0.05, [d["baseline"][k] for k, _ in metrics_keys], s=14, color="#1E293B", zorder=5)
            ax_sim.scatter(xs + bw / 2 + (rep_i - 1) * 0.05, [d["candidate"][k] for k, _ in metrics_keys], s=14, color="#0369A1", zorder=5)
        ax_sim.set_xticks(xs)
        ax_sim.set_xticklabels([lbl for _, lbl in metrics_keys], fontsize=8.2)
        ax_sim.set_ylim(0, 1.15)
        ax_sim.set_ylabel("Score vs NanoSim Ground Truth")
        ax_sim.set_title(f"D. SIRV Ground-Truth Accuracy Across simulated1..{len(sim_docs)} (+ GENCODE v46)", fontweight="bold")
        ax_sim.grid(axis="y", linestyle="--", alpha=0.35)
        ax_sim.legend(loc="lower right", frameon=True, fontsize=8.0)
        for xi, bv, cv in zip(xs, b_vals, c_vals):
            ax_sim.text(xi - bw / 2, bv + 0.015, f"{bv:.4f}", ha="center", va="bottom", fontsize=7.0, rotation=45, color="#334155")
            ax_sim.text(xi + bw / 2, cv + 0.015, f"{cv:.4f}", ha="center", va="bottom", fontsize=7.0, rotation=45, color="#0369A1", fontweight="bold")

    fig.suptitle(
        "Quantitative Correlation, Read/Model Parity & NanoSim Ground-Truth Accuracy Across Full Cohort",
        fontsize=12.8,
        fontweight="bold",
        y=0.985,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_fig4_pss_profiles(bench_root: Path, out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.0, "axes.edgecolor": "#CBD5E1"})
    fig, axes = plt.subplots(2, 2, figsize=(15.8, 9.6))
    panels = [
        (axes[0, 0], "g1_encff400bqq_t8_abba3", "rep2", "A. Tier 1 ENCFF400BQQ (-t 8) Cold vs Warm .iqidx", "g1_encff400bqq_t8_warm_abba3"),
        (axes[0, 1], "g3_shrbmx_rep2_v46_t32", "rep1", "B. Tier 3 shRBMX_rep2 (20.44M ONT + v46, -t 32)", None),
        (axes[1, 0], "g1_cdna_pb_manatee_trichechus_manatus_encff589ldn_t16", "rep1", "C. Manatee ENCFF589LDN De Novo (15,608 Contigs, -92% PSS)", None),
        (axes[1, 1], "g3_encff089amf_t32", "rep1", "D. Tier 3 ENCFF089AMF (16.42M Capped ONT, -t 32)", None),
    ]
    for ax, slug, rep, title, warm_slug in panels:
        bo = bench_root / slug / f"baseline.{rep}.pss.tsv"
        bn = bench_root / slug / f"candidate.{rep}.pss.tsv"
        if bo.exists() and bn.exists():
            t_o, r_o, p_o = _load_pss_trace(bo)
            t_n, r_n, p_n = _load_pss_trace(bn)
            ax.plot(t_o, r_o, color=C_ORIG, linewidth=1.3, linestyle="--", alpha=0.68, label=f"Orig RSS (peak {r_o.max():.2f} GiB)")
            ax.plot(t_o, p_o, color=C_ORIG, linewidth=2.0, label=f"Orig PSS ({t_o[-1]:.1f}m, peak {p_o.max():.2f} GiB)")
            ax.plot(t_n, r_n, color=C_NEXT, linewidth=1.3, linestyle="--", alpha=0.80, label=f"Next RSS (peak {r_n.max():.2f} GiB)")
            ax.plot(t_n, p_n, color=C_ACCENT, linewidth=2.2, label=f"Next PSS ({t_n[-1]:.1f}m, peak {p_n.max():.2f} GiB)")
            if warm_slug and (bench_root / warm_slug / f"candidate.{rep}.pss.tsv").exists():
                t_w, _, p_w = _load_pss_trace(bench_root / warm_slug / f"candidate.{rep}.pss.tsv")
                ax.plot(t_w, p_w, color=C_INDIGO, linewidth=2.0, label=f"Warm .iqidx PSS ({t_w[-1]:.2f}m, peak {p_w.max():.2f} GiB)")
            ax.set_ylim(-0.5, max(r_o.max(), p_o.max(), r_n.max()) * 1.15)
        ax.set_xlabel("Elapsed Wall-Clock Time (minutes)")
        ax.set_ylabel("Process-Tree Memory (GiB)")
        ax.set_title(title, fontweight="bold")
        ax.grid(True, linestyle="--", alpha=0.35)
        ax.legend(loc="upper right", frameon=True, framealpha=0.95, fontsize=7.8)

    fig.suptitle("Live Process-Tree Memory Trajectories Over Execution Time (/proc/<pid>/smaps_rollup)", fontsize=13, fontweight="bold", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_fig5_species_cellline_heatmap(cohort: dict, out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.0, "axes.edgecolor": "#CBD5E1"})
    paired = [r for r in cohort["rows"] if r.get("is_paired_engine_cmp")]
    row_protos = ["cDNA_PacBio", "capped_PacBio", "directRNA", "cDNA_ONT", "capped_ONT", "R2C2_ONT", "R2C2_ONT_rmdup", "RBMX_KD"]
    col_cells = [
        ("Human", "H1_DE_cell", "Human\nH1_DE"),
        ("Human", "H1_ES_cell", "Human\nH1_ES"),
        ("Human", "WTC11_cell", "Human\nWTC11"),
        ("Human", "H1_DE_cell_seqkitrmdup", "Human\nH1_DE\nrmdup"),
        ("Human", "H1_ES_cell_seqkitrmdup", "Human\nH1_ES\nrmdup"),
        ("Human", "WTC11_RBMX_KD", "Human\nRBMX KD"),
        ("Mouse", "F121_cell", "Mouse\nF121"),
        ("Manatee", "Trichechus_manatus", "Manatee\n(De Novo)"),
    ]
    # Keep only cells that appear in paired rows
    present_cells = [(sp, cl, lbl) for sp, cl, lbl in col_cells if any(r["species"] == sp and r["cell_line"] == cl for r in paired)]
    present_protos = [p for p in row_protos if any(r["protocol"] == p for r in paired)]

    spd_mat = np.full((len(present_protos), len(present_cells)), np.nan)
    pss_mat = np.full((len(present_protos), len(present_cells)), np.nan)
    spd_txt = [["" for _ in present_cells] for _ in present_protos]
    cor_txt = [["" for _ in present_cells] for _ in present_protos]

    for i, p in enumerate(present_protos):
        for j, (sp, cl, _) in enumerate(present_cells):
            items = [r for r in paired if r["protocol"] == p and r["species"] == sp and r["cell_line"] == cl]
            if not items:
                continue
            sm = float(np.median([x["speedup_x"] for x in items]))
            pm = float(np.median([x["pss_delta_pct"] for x in items if x["pss_delta_pct"] is not None]))
            rs = [x["primary_tx_or_disc_r"] for x in items if x.get("primary_tx_or_disc_r") is not None]
            rmin = min(rs) if rs else 1.0
            spd_mat[i, j] = sm
            pss_mat[i, j] = pm
            spd_txt[i][j] = f"{sm:.1f}x\n(n={len(items)})"
            cor_txt[i][j] = f"PSS -{pm:.0f}%\nr>={rmin:.5f}"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.0, 6.8), sharey=True)
    cmap1 = plt.cm.YlGnBu.copy()
    cmap1.set_bad("#F1F5F9")
    im1 = ax1.imshow(spd_mat, cmap=cmap1, aspect="auto", vmin=1.0, vmax=max(20.0, float(np.nanmax(spd_mat))))
    ax1.set_xticks(np.arange(len(present_cells)))
    ax1.set_xticklabels([c[2] for c in present_cells], fontsize=8.3)
    ax1.set_yticks(np.arange(len(present_protos)))
    ax1.set_yticklabels(present_protos, fontsize=9.0, fontweight="medium")
    ax1.set_title("A. Median Wall-Clock Speedup (x) by Protocol x Species / Cell Line", fontweight="bold")
    for i in range(len(present_protos)):
        for j in range(len(present_cells)):
            if spd_txt[i][j]:
                col = "white" if spd_mat[i, j] > 12 else "#0F172A"
                ax1.text(j, i, spd_txt[i][j], ha="center", va="center", fontsize=8.0, fontweight="bold", color=col)
            else:
                ax1.text(j, i, "—", ha="center", va="center", fontsize=9, color="#94A3B8")
    fig.colorbar(im1, ax=ax1, fraction=0.035, pad=0.02, label="Median Speedup (x)")

    cmap2 = plt.cm.Greens.copy()
    cmap2.set_bad("#F1F5F9")
    im2 = ax2.imshow(pss_mat, cmap=cmap2, aspect="auto", vmin=0.0, vmax=95.0)
    ax2.set_xticks(np.arange(len(present_cells)))
    ax2.set_xticklabels([c[2] for c in present_cells], fontsize=8.3)
    ax2.set_title("B. Median Peak PSS Savings (%) & Minimum Correlation r", fontweight="bold")
    for i in range(len(present_protos)):
        for j in range(len(present_cells)):
            if cor_txt[i][j]:
                col = "white" if pss_mat[i, j] > 58 else "#0F172A"
                ax2.text(j, i, cor_txt[i][j], ha="center", va="center", fontsize=7.6, fontweight="bold", color=col)
            else:
                ax2.text(j, i, "—", ha="center", va="center", fontsize=9, color="#94A3B8")
    fig.colorbar(im2, ax=ax2, fraction=0.035, pad=0.02, label="Peak PSS Savings (%)")

    fig.suptitle("Protocol x Species x Cell-Line Matrix: Speedup, Memory Compression & Expression Parity", fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.subplots_adjust(wspace=0.22)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_fig6_features_and_vis(cohort: dict, bench_root: Path, out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.0, "axes.edgecolor": "#CBD5E1"})
    fig, axes = plt.subplots(2, 2, figsize=(15.8, 9.6))
    ax_cage, ax_res, ax_misc, ax_vis = axes.ravel()

    # Panel A: CAGE TSS discovery across RBMX_KD
    cage_rows = [r for r in cohort["rows"] if r["protocol"] == "RBMX_KD_CAGE" and r.get("cage_info")]
    if cage_rows:
        xc = np.arange(len(cage_rows))
        tss = [r["cage_info"]["tss_prediction_lines"] / 1000.0 for r in cage_rows]
        mods = [r["cage_info"]["transcript_models_count"] / 1000.0 for r in cage_rows]
        walls = [r["next_wall_s"] / 60.0 for r in cage_rows]
        bw = 0.36
        ax_cage.bar(xc - bw / 2, tss, width=bw, color=C_NEXT, label="TSS Predictions (x1,000 lines)")
        ax_cage.bar(xc + bw / 2, mods, width=bw, color=C_ACCENT, label="Transcript Models (x1,000)")
        ax_cage.set_xticks(xc)
        ax_cage.set_xticklabels([f"{r['sample_id']}\n({w:.1f}m)" for r, w in zip(cage_rows, walls)])
        ax_cage.set_ylabel("Count (x1,000)")
        ax_cage.set_ylim(0, max(tss) * 1.22)
        ax_cage.set_title("A. CAGE-Assisted TSS Discovery (isoquant-next --cage + v46)", fontweight="bold")
        ax_cage.grid(axis="y", linestyle="--", alpha=0.35)
        ax_cage.legend(loc="upper right", frameon=True, fontsize=8.0)
        for xi, tv, mv in zip(xc, tss, mods):
            ax_cage.text(xi - bw / 2, tv + 0.3, f"{tv:.1f}k", ha="center", va="bottom", fontsize=8.0, fontweight="bold", color="#0369A1")
            ax_cage.text(xi + bw / 2, mv + 0.3, f"{mv:.1f}k", ha="center", va="bottom", fontsize=8.0, fontweight="bold", color="#065F46")

    # Panel B: Resume & Illumina & Multi-Sample YAML Parity
    special_labels = {
        "g1_multisample_yaml_t16": "MultiSample YAML\n(H1_ES 2-rep)",
        "g1_cross_cell_yaml_t16": "MultiSample YAML\n(Cross-Cell)",
        "g2_rbmx_yaml_t32": "MultiSample YAML\n(RBMX 3-cond)",
        "g4_chr9_illumina_t8": "Illumina Hybrid\n(Sim chr9)",
        "g4_resume_visualize_t16": "Resume & Vis\n(ENCFF400BQQ)",
    }
    special = []
    for r in cohort["rows"]:
        if r["protocol"] in ("Resume_Visualize", "Illumina_Hybrid", "MultiSample_YAML"):
            lbl = special_labels.get(r["slug"], f"{r['protocol'].replace('_', chr(10))}\n({r['sample_id'][:12]})")
            special.append((lbl, r))
    if special:
        xs = np.arange(len(special))
        vals_r = [r.get("primary_tx_or_disc_r") or 1.0 for _, r in special]
        vals_j = [r.get("transcript_models_jaccard") or 1.0 for _, r in special]
        bw = 0.35
        ax_res.bar(xs - bw / 2, vals_r, width=bw, color=C_INDIGO, label="Tx Pearson r / Exact Match")
        ax_res.bar(xs + bw / 2, vals_j, width=bw, color=C_WARN, label="Transcript Models Jaccard")
        ax_res.set_xticks(xs)
        ax_res.set_xticklabels([l for l, _ in special], fontsize=7.6)
        ax_res.set_ylim(0.985, 1.003)
        ax_res.set_ylabel("Parity Score")
        ax_res.set_title("B. --resume, --illumina_bam & --yaml Multi-Sample Parity", fontweight="bold")
        ax_res.grid(axis="y", linestyle="--", alpha=0.35)
        ax_res.legend(loc="lower right", frameon=True, fontsize=8.0)
        for xi, vr, vj, (_, r) in zip(xs, vals_r, vals_j, special):
            ax_res.text(xi - bw / 2, vr + 0.0005, f"{vr:.5f}", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=C_INDIGO)
            ax_res.text(xi + bw / 2, vj + 0.0005, f"{vj:.4f}", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color="#B45309")

    # Panel C: Sequin & FASTQ & Reference Version Contrasts
    contrasts = [
        r for r in cohort["rows"]
        if r["protocol"] in ("Sequin_SIRV4", "FASTQ_directRNA", "cDNA_PacBio_v46", "cDNA_PacBio_vM36")
    ]
    if contrasts:
        xm = np.arange(len(contrasts))
        spds = [r["speedup_x"] or 0.0 for r in contrasts]
        pss_s = [r["pss_delta_pct"] or 0.0 for r in contrasts]
        ax_misc.bar(xm, spds, width=0.52, color=C_NEXT, label="Wall-Clock Speedup (x)")
        ax_misc.set_xticks(xm)
        ax_misc.set_xticklabels([f"{r['protocol']}\n({r['sample_id']})" for r in contrasts], fontsize=7.8)
        ax_misc.set_ylabel("Speedup Factor (x)")
        ax_misc.set_ylim(0, max(spds) * 1.28)
        ax_misc.set_title("C. Reference Contrasts (v46/vM36), Sequin & FASTQ De Novo", fontweight="bold")
        ax_misc.grid(axis="y", linestyle="--", alpha=0.35)
        for xi, sp, ps, r in zip(xm, spds, pss_s, contrasts):
            pr = r.get("primary_tx_or_disc_r") or 1.0
            ps_str = f"-{ps:.1f}%" if ps >= 0 else f"+{abs(ps):.1f}%"
            ax_misc.text(xi, sp + max(spds) * 0.02, f"{sp:.2f}x\nPSS {ps_str}\nr={pr:.6f}", ha="center", va="bottom", fontsize=7.4, fontweight="bold", color="#0F172A")

    # Panel D: Native isoquant-next-visualize sample plot
    vis_dir = bench_root / "g4_resume_visualize_t16" / "vis_png"
    pngs = sorted(vis_dir.glob("*.png")) if vis_dir.exists() else []
    if pngs:
        img = mpimg.imread(str(pngs[0]))
        ax_vis.imshow(img)
        ax_vis.axis("off")
        ax_vis.set_title(f"D. Native isoquant-next-visualize Output ({pngs[0].stem})", fontweight="bold")
    else:
        ax_vis.axis("off")
        ax_vis.set_title("D. isoquant-next-visualize (pending)", fontweight="bold")

    fig.suptitle("Gate G2–G4 Feature Closures: CAGE, --resume, --illumina_bam, --yaml, Reference Contrasts & Visualization", fontsize=12.6, fontweight="bold", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot full-cohort figures from cohort_all_results.json.")
    ap.add_argument("--cohort-json", type=Path, required=True)
    ap.add_argument("--bench-root", type=Path, required=True)
    ap.add_argument("--scatter-dir", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    cohort = json.loads(args.cohort_json.read_text(encoding="utf-8"))

    plot_fig1_speedup(cohort, args.out_dir / "fig1_speedup_and_walltime.png")
    plot_fig2_memory_and_scaling(cohort, args.out_dir / "fig2_memory_and_scaling.png")
    plot_fig3_correlation_and_sirv(cohort, args.scatter_dir, args.bench_root, args.out_dir / "fig3_correlation_and_accuracy.png")
    plot_fig4_pss_profiles(args.bench_root, args.out_dir / "fig4_pss_memory_profiles.png")
    plot_fig5_species_cellline_heatmap(cohort, args.out_dir / "fig5_species_cellline_heatmap.png")
    plot_fig6_features_and_vis(cohort, args.bench_root, args.out_dir / "fig6_features_cage_resume_vis.png")
    print(f"Generated 6 cohort figures in {args.out_dir}")


if __name__ == "__main__":
    main()
