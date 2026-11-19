#!/usr/bin/env python3
"""Generate publication-quality Phase 5 performance, memory, and correlation charts."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np


BENCHMARKS_ORDER: List[Tuple[str, str, str]] = [
    ("g1_encff400bqq_t1", "G1: ENCFF400BQQ (-t 1)", "G1"),
    ("g1_encff400bqq_t8_abba3", "G1: ENCFF400BQQ (3x ABBA, -t 8)", "G1"),
    ("g1_encff400bqq_t16", "G1: ENCFF400BQQ (-t 16)", "G1"),
    ("g1_encff400bqq_t8_warm_abba3", "G1: ENCFF400BQQ Warm (3x ABBA, -t 8)", "G1"),
    ("g1_encff684yoo_t16", "G1: ENCFF684YOO CapPacBio (-t 16)", "G1"),
    ("g1_encff339fmq_t16", "G1: ENCFF339FMQ CapPacBio (-t 16)", "G1"),
    ("g1_multisample_yaml_t16", "G1: Multi-Sample YAML (-t 16)", "G1"),
    ("g2_encff436gkz_t16", "G2: ENCFF436GKZ cDNA (-t 16)", "G2"),
    ("g2_encff846rid_drna_t16", "G2: ENCFF846RID dRNA (-t 16)", "G2"),
    ("g2_encff207kvo_drna_t16", "G2: ENCFF207KVO dRNA (-t 16)", "G2"),
    ("g2_con_rep1_v46_t16", "G2: con_rep1 6.7M ONT+v46 (-t 16)", "G2"),
    ("g3_simulated1_t16", "G3: simulated1 8.8M+v46 (-t 16)", "G3"),
    ("g3_encff598yqo_mouse_t16", "G3: ENCFF598YQO 10.1M Mouse (-t 16)", "G3"),
    ("g3_encff089amf_t32", "G3: ENCFF089AMF 16.4M ONT (-t 32)", "G3"),
    ("g3_shrbmx_rep2_v46_t32", "G3: shRBMX_rep2 20.4M+v46 (-t 32)", "G3"),
    ("g4_encff008hmy_fastq_t16", "G4: ENCFF008HMY.fastq -> minimap2 (-t 16)", "G4"),
]

C_ORIG = "#64748B"
C_NEXT = "#0284C7"
C_ACCENT = "#059669"
C_WARN = "#D97706"
C_RSS_ORIG = "#94A3B8"
C_RSS_NEXT = "#38BDF8"


def _fmt_duration(sec: float) -> str:
    if sec >= 3600:
        return f"{sec / 3600:.2f}h"
    if sec >= 120:
        return f"{sec / 60:.1f}m"
    return f"{sec:.1f}s"


def _load_all(bench_root: Path) -> List[dict]:
    rows = []
    for slug, label, tier in BENCHMARKS_ORDER:
        d = bench_root / slug
        bs = json.loads((d / "benchmark_summary.json").read_text())["summary"]
        base = bs["baseline"]
        cand = bs["candidate"]
        acc_files = sorted(d.glob("accuracy*.json"))
        accs = [json.loads(p.read_text()) for p in acc_files]
        rows.append(
            {
                "slug": slug,
                "label": label,
                "tier": tier,
                "orig_wall": float(base["wall_median_s"]),
                "next_wall": float(cand["wall_median_s"]),
                "speedup": float(base["wall_median_s"]) / float(cand["wall_median_s"]),
                "orig_pss_gib": float(base["pss_median_kib"]) / 1048576.0,
                "next_pss_gib": float(cand["pss_median_kib"]) / 1048576.0,
                "orig_rss_gib": float(base["rss_median_kib"]) / 1048576.0,
                "next_rss_gib": float(cand["rss_median_kib"]) / 1048576.0,
                "accs": accs,
            }
        )
    return rows


def plot_fig1_speedup(rows: List[dict], out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": "#CBD5E1"})
    n = len(rows)
    y = np.arange(n)[::-1]
    h = 0.36

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(15.5, 8.2), gridspec_kw={"width_ratios": [1.35, 1.0]}, sharey=True
    )

    orig_w = [r["orig_wall"] for r in rows]
    next_w = [r["next_wall"] for r in rows]
    spd = [r["speedup"] for r in rows]
    labels = [r["label"] for r in rows]

    ax1.barh(y + h / 2, orig_w, height=h, color=C_ORIG, label="Original isoquant.py (Baseline)")
    ax1.barh(y - h / 2, next_w, height=h, color=C_NEXT, label="Refactored isoquant-next (Scheme B)")
    ax1.set_xscale("log")
    ax1.set_xlim(20, 48000)
    ax1.set_ylim(-0.7, n - 0.3)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontweight="medium")
    ax1.set_xlabel("Wall-Clock Runtime (seconds, log scale)")
    ax1.set_title("A. End-to-End Runtime (Baseline vs Refactored)", fontweight="bold", pad=10)
    ax1.grid(axis="x", linestyle="--", alpha=0.4)
    handles, lbls = ax1.get_legend_handles_labels()
    fig.legend(handles, lbls, loc="upper center", bbox_to_anchor=(0.56, 0.945), ncol=2, frameon=True, framealpha=0.95)

    for yi, ow, nw in zip(y, orig_w, next_w):
        ax1.text(ow * 1.12, yi + h / 2, _fmt_duration(ow), va="center", ha="left", fontsize=7.8, color="#334155")
        ax1.text(nw * 1.12, yi - h / 2, _fmt_duration(nw), va="center", ha="left", fontsize=7.8, color="#0369A1", fontweight="bold")

    bar_colors = [C_ACCENT if s >= 10.0 else C_NEXT for s in spd]
    ax2.barh(y, spd, height=0.58, color=bar_colors, alpha=0.92)
    ax2.axvline(1.0, color="#94A3B8", linestyle="--", linewidth=1.1, label="1.0x Parity")
    ax2.set_xlim(0, max(spd) * 1.18)
    ax2.set_xlabel("Speedup Factor (x faster than isoquant.py)")
    ax2.set_title("B. Measured Speedup Across Tiers 1–4 (2.93x – 31.02x)", fontweight="bold", pad=10)
    ax2.grid(axis="x", linestyle="--", alpha=0.4)
    for yi, s in zip(y, spd):
        ax2.text(s + 0.45, yi, f"{s:.2f}x", va="center", ha="left", fontsize=8.4, fontweight="bold", color="#0F172A")

    for sep in [8.5, 4.5, 0.5]:
        ax1.axhline(sep, color="#E2E8F0", linewidth=0.9, linestyle=":")
        ax2.axhline(sep, color="#E2E8F0", linewidth=0.9, linestyle=":")

    fig.suptitle(
        "IsoQuant Scheme B Refactoring: End-to-End Speedup Across 16 Paired Slurm Benchmarks",
        fontsize=13,
        fontweight="bold",
        y=0.985,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_fig2_memory_and_scaling(rows: List[dict], out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.2, "axes.edgecolor": "#CBD5E1"})
    fig = plt.figure(figsize=(15.5, 9.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.35, 1.0])
    ax_pss = fig.add_subplot(gs[0, :])
    ax_rss = fig.add_subplot(gs[1, 0])
    ax_scale = fig.add_subplot(gs[1, 1])

    x = np.arange(len(rows))
    w = 0.38
    short_names = [
        "400BQQ\n-t1",
        "400BQQ\nABBA -t8",
        "400BQQ\n-t16",
        "400BQQ\nWarm -t8",
        "684YOO\n-t16",
        "339FMQ\n-t16",
        "Multi\nYAML",
        "436GKZ\n-t16",
        "846RID\ndRNA",
        "207KVO\ndRNA",
        "con_rep1\nv46",
        "sim1\nv46",
        "598YQO\nMouse",
        "089AMF\n16M -t32",
        "shRBMX\n20M -t32",
        "008HMY\nFASTQ*",
    ]
    opss = [r["orig_pss_gib"] for r in rows]
    npss = [r["next_pss_gib"] for r in rows]

    ax_pss.bar(x - w / 2, opss, width=w, color=C_ORIG, label="Original isoquant.py Peak Tree PSS (GiB)")
    ax_pss.bar(x + w / 2, npss, width=w, color=C_NEXT, label="Refactored isoquant-next Peak Tree PSS (GiB)")
    ax_pss.set_xticks(x)
    ax_pss.set_xticklabels(short_names, fontsize=8.0)
    ax_pss.set_ylabel("Peak Process-Tree PSS (GiB)")
    ax_pss.set_ylim(0, max(opss) * 1.18)
    ax_pss.set_title(
        "A. Peak Process-Tree Memory (PSS via /proc/<pid>/smaps_rollup) Across 16 Paired Runs (*G4 FASTQ dominated by minimap2 ~21.6 GiB genome index)",
        fontweight="bold",
    )
    ax_pss.grid(axis="y", linestyle="--", alpha=0.4)
    ax_pss.legend(loc="upper left", frameon=True)

    for xi, op, np_ in zip(x, opss, npss):
        pct = (op - np_) / op * 100.0
        txt = f"-{pct:.0f}%" if pct >= 1.0 else "~0%"
        col = C_ACCENT if pct >= 50.0 else ("#0369A1" if pct > 5.0 else "#64748B")
        ax_pss.text(xi + w / 2, np_ + 0.4, txt, ha="center", va="bottom", fontsize=7.7, fontweight="bold", color=col)

    heavy_slugs = [
        ("g1_encff400bqq_t16", "400BQQ\n-t16"),
        ("g2_con_rep1_v46_t16", "con_rep1\nv46 -t16"),
        ("g3_simulated1_t16", "simulated1\n8.8M -t16"),
        ("g3_encff089amf_t32", "ENCFF089AMF\n16.4M -t32"),
        ("g3_shrbmx_rep2_v46_t32", "shRBMX_rep2\n20.4M -t32"),
    ]
    by_slug = {r["slug"]: r for r in rows}
    xh = np.arange(len(heavy_slugs))
    orss = [by_slug[s]["orig_rss_gib"] for s, _ in heavy_slugs]
    nrss = [by_slug[s]["next_rss_gib"] for s, _ in heavy_slugs]
    ax_rss.bar(xh - w / 2, orss, width=w, color=C_ORIG, label="Original Peak Tree RSS (GiB)")
    ax_rss.bar(xh + w / 2, nrss, width=w, color=C_NEXT, label="Refactored Peak Tree RSS (GiB)")
    ax_rss.set_xticks(xh)
    ax_rss.set_xticklabels([lbl for _, lbl in heavy_slugs], fontsize=8.5)
    ax_rss.set_ylabel("Peak Tree RSS (GiB)")
    ax_rss.set_ylim(0, max(orss) * 1.18)
    ax_rss.set_title("B. High-Concurrency RSS Savings (Copy-on-Write vs Shared mmap .iqidx)", fontweight="bold")
    ax_rss.grid(axis="y", linestyle="--", alpha=0.4)
    ax_rss.legend(loc="upper left", frameon=True, fontsize=8.2)
    for xi, orv, nrv in zip(xh, orss, nrss):
        red = (orv - nrv) / orv * 100.0
        ax_rss.text(xi - w / 2, orv + 0.8, f"{orv:.1f}G", ha="center", va="bottom", fontsize=7.8, color="#334155")
        ax_rss.text(xi + w / 2, nrv + 0.8, f"{nrv:.1f}G\n(-{red:.0f}%)", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color=C_ACCENT)

    threads = [1, 8, 16]
    scale_slugs = ["g1_encff400bqq_t1", "g1_encff400bqq_t8_abba3", "g1_encff400bqq_t16"]
    o_times = [by_slug[s]["orig_wall"] / 60.0 for s in scale_slugs]
    n_times = [by_slug[s]["next_wall"] / 60.0 for s in scale_slugs]
    warm_min = by_slug["g1_encff400bqq_t8_warm_abba3"]["next_wall"] / 60.0

    ax_scale.plot(threads, o_times, "o-", color=C_ORIG, linewidth=2.2, markersize=7, label="Original isoquant.py (Cold GTF -> DB)")
    ax_scale.plot(threads, n_times, "s-", color=C_NEXT, linewidth=2.2, markersize=7, label="Refactored isoquant-next (Cold GTF -> .iqidx)")
    ax_scale.scatter([8], [warm_min], marker="*", s=200, color=C_ACCENT, zorder=5, label=f"isoquant-next Warm .iqidx (-t 8: {warm_min:.2f}m, 31.0x)")
    ax_scale.set_yscale("log")
    ax_scale.set_ylim(0.55, 280)
    ax_scale.set_xlim(-1.0, 17.8)
    ax_scale.set_xticks(threads)
    ax_scale.set_xticklabels(["-t 1", "-t 8", "-t 16"])
    ax_scale.set_ylabel("Wall-Clock Runtime (minutes, log scale)")
    ax_scale.set_title("C. Thread Scaling on ENCFF400BQQ (-t 1 -> 8 -> 16 + Warm .iqidx)", fontweight="bold")
    ax_scale.grid(True, which="both", linestyle="--", alpha=0.35)
    ax_scale.legend(loc="upper right", frameon=True, fontsize=8.0)

    for th, ot, nt in zip(threads, o_times, n_times):
        ha = "left" if th == 1 else "center"
        x_t = th - 0.1 if th == 1 else th
        ax_scale.text(x_t, ot * 1.18, f"{ot:.1f}m", ha=ha, va="bottom", fontsize=8.2, color="#334155")
        ax_scale.text(x_t, nt * 1.24, f"{nt:.2f}m ({ot/nt:.1f}x)", ha=ha, va="bottom", fontsize=8.0, fontweight="bold", color="#0369A1")
    ax_scale.text(8, warm_min * 0.72, f"Warm .iqidx: {warm_min:.2f}m (31.0x)", ha="center", va="top", fontsize=8.0, fontweight="bold", color=C_ACCENT)

    fig.suptitle("Process-Tree Memory Footprint (PSS / RSS) & Parallel Thread Scaling", fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


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


def plot_fig3_correlation(rows: List[dict], scatter_dir: Path, bench_root: Path, out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.2, "axes.edgecolor": "#CBD5E1"})
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 11.2))
    ax_g, ax_t, ax_par, ax_sim = axes.ravel()

    # Panel A & B: Gene and Transcript correlation scatter on ENCFF400BQQ chr1
    og = _load_counts(scatter_dir / "orig_gene.tsv")
    ng = _load_counts(scatter_dir / "next_gene.tsv")
    g_keys = sorted(set(og) | set(ng))
    gx = np.array([og.get(k, 0.0) for k in g_keys])
    gy = np.array([ng.get(k, 0.0) for k in g_keys])
    mask_g = (gx > 0) | (gy > 0)
    g_r = float(np.corrcoef(gx, gy)[0, 1])
    g_exact = float(np.mean(np.abs(gx - gy) <= 1e-6))

    lx, ly = np.log10(gx[mask_g] + 1.0), np.log10(gy[mask_g] + 1.0)
    lim_g = max(lx.max(), ly.max()) * 1.05
    ax_g.scatter(lx, ly, s=16, alpha=0.55, color=C_NEXT, edgecolors="none", label=f"Expressed genes (n={mask_g.sum():,})")
    ax_g.plot([0, lim_g], [0, lim_g], "--", color="#0F172A", linewidth=1.1, label="Identity y = x")
    ax_g.set_xlim(-0.08, lim_g)
    ax_g.set_ylim(-0.08, lim_g)
    ax_g.set_xlabel("Original isoquant.py: log10(Gene Count + 1)")
    ax_g.set_ylabel("Refactored isoquant-next: log10(Gene Count + 1)")
    ax_g.set_title("A. Per-Gene Quantification Correlation (ENCFF400BQQ chr21+22)", fontweight="bold")
    ax_g.grid(True, linestyle="--", alpha=0.35)
    ax_g.text(
        0.04, 0.88,
        f"Pearson r (raw) = {g_r:.6f}\nExact match (all {len(g_keys):,}) = {g_exact*100:.2f}%\nMax |diff| = {np.max(np.abs(gx - gy)):.4g}",
        transform=ax_g.transAxes, va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8FAFC", edgecolor="#CBD5E1"),
        fontsize=8.6,
    )
    ax_g.legend(loc="lower right", frameon=True)

    ot = _load_counts(scatter_dir / "orig_tx.tsv")
    nt = _load_counts(scatter_dir / "next_tx.tsv")
    t_keys = sorted(set(ot) | set(nt))
    tx = np.array([ot.get(k, 0.0) for k in t_keys])
    ty = np.array([nt.get(k, 0.0) for k in t_keys])
    mask_t = (tx > 0) | (ty > 0)
    t_r = float(np.corrcoef(tx, ty)[0, 1])
    t_exact = float(np.mean(np.abs(tx - ty) <= 1e-6))

    ltx, lty = np.log10(tx[mask_t] + 1.0), np.log10(ty[mask_t] + 1.0)
    lim_t = max(ltx.max(), lty.max()) * 1.05
    ax_t.scatter(ltx, lty, s=14, alpha=0.48, color=C_ACCENT, edgecolors="none", label=f"Expressed transcripts (n={mask_t.sum():,})")
    ax_t.plot([0, lim_t], [0, lim_t], "--", color="#0F172A", linewidth=1.1, label="Identity y = x")
    ax_t.set_xlim(-0.08, lim_t)
    ax_t.set_ylim(-0.08, lim_t)
    ax_t.set_xlabel("Original isoquant.py: log10(Transcript Count + 1)")
    ax_t.set_ylabel("Refactored isoquant-next: log10(Transcript Count + 1)")
    ax_t.set_title("B. Per-Transcript Quantification Correlation (ENCFF400BQQ chr21+22)", fontweight="bold")
    ax_t.grid(True, linestyle="--", alpha=0.35)
    ax_t.text(
        0.04, 0.88,
        f"Pearson r (raw) = {t_r:.6f}\nExact match (all {len(t_keys):,}) = {t_exact*100:.2f}%\nMax |diff| = {np.max(np.abs(tx - ty)):.4g}",
        transform=ax_t.transAxes, va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8FAFC", edgecolor="#CBD5E1"),
        fontsize=8.6,
    )
    ax_t.legend(loc="lower right", frameon=True)

    # Panel C: Quantitative Parity Metrics Across All 16 Paired Runs
    par_labels, tx_rs, ri_agrees, ext_js, tm_js = [], [], [], [], []
    for r in rows:
        for a in r["accs"]:
            suf = f" ({a['prefix']})" if len(r["accs"]) > 1 else ""
            par_labels.append(r["label"].split(": ", 1)[1] + suf)
            tx_rs.append(a["transcript_counts"]["pearson_r"])
            ri_agrees.append(a["read_info"]["full_assignment_agreement"])
            ext_js.append(a["extended_annotation_gtf"]["jaccard"])
            tm_js.append(a["transcript_models_gtf"]["jaccard"])

    yp = np.arange(len(par_labels))[::-1]
    ax_par.scatter(tx_rs, yp, marker="o", s=48, color="#1D4ED8", label="Transcript Counts Pearson r (min 0.99999985)", zorder=4)
    ax_par.scatter(ri_agrees, yp, marker="s", s=44, color=C_ACCENT, label="Read Assignment Full Agreement (min 0.999606)", zorder=3)
    ax_par.scatter(ext_js, yp, marker="D", s=38, color="#7C3AED", label="Extended Annotation GTF Jaccard (min 0.999922)", zorder=3)
    ax_par.scatter(tm_js, yp, marker="^", s=44, color=C_WARN, label="Transcript Models GTF Jaccard (min 0.996655)", zorder=2)
    ax_par.axvline(1.0, color="#0F172A", linestyle=":", linewidth=1.0)
    ax_par.set_xlim(0.9932, 1.00045)
    ax_par.set_yticks(yp)
    ax_par.set_yticklabels(par_labels, fontsize=7.8)
    ax_par.set_xlabel("Correlation / Agreement / Structural Jaccard Score (Zoomed 0.9932 – 1.0000)")
    ax_par.set_title("C. Output Correlation & Structural Parity Across All 17 Evaluation Pairs", fontweight="bold")
    ax_par.grid(True, linestyle="--", alpha=0.35)
    ax_par.legend(loc="upper left", frameon=True, framealpha=0.95, fontsize=7.6)

    # Panel D: Gate G3 simulated1 Ground-Truth Recall/Precision
    sim_eval = json.loads((bench_root / "g3_simulated1_t16" / "eval_simulated1.json").read_text())
    b_ev, c_ev = sim_eval["baseline"], sim_eval["candidate"]
    metrics_keys = [
        ("read_precision", "Read\nPrecision"),
        ("read_recall", "Read\nRecall"),
        ("model_known_over_all_precision", "Model Known\nPrecision"),
        ("model_recall_ge2_reads", "Model Recall\n(>=2 reads)"),
        ("quant_precision", "Quant\nPrecision"),
        ("quant_recall_ge2_reads", "Quant Recall\n(>=2 reads)"),
    ]
    xs = np.arange(len(metrics_keys))
    bw = 0.36
    b_vals = [b_ev[k] for k, _ in metrics_keys]
    c_vals = [c_ev[k] for k, _ in metrics_keys]
    ax_sim.bar(xs - bw / 2, b_vals, width=bw, color=C_ORIG, label="Original isoquant.py")
    ax_sim.bar(xs + bw / 2, c_vals, width=bw, color=C_NEXT, label="Refactored isoquant-next")
    ax_sim.set_xticks(xs)
    ax_sim.set_xticklabels([lbl for _, lbl in metrics_keys], fontsize=8.3)
    ax_sim.set_ylim(0, 1.12)
    ax_sim.set_ylabel("Score vs NanoSim Ground Truth")
    ax_sim.set_title("D. Ground-Truth Accuracy on simulated1 (8.77M NanoSim Reads + v46)", fontweight="bold")
    ax_sim.grid(axis="y", linestyle="--", alpha=0.35)
    ax_sim.legend(loc="lower right", frameon=True)
    for xi, bv, cv in zip(xs, b_vals, c_vals):
        ax_sim.text(xi - bw / 2, bv + 0.015, f"{bv:.4f}", ha="center", va="bottom", fontsize=7.0, rotation=45, color="#334155")
        ax_sim.text(xi + bw / 2, cv + 0.015, f"{cv:.4f}", ha="center", va="bottom", fontsize=7.0, rotation=45, color="#0369A1", fontweight="bold")

    fig.suptitle(
        "Quantitative Output Correlation, Read-Level Agreement, and Ground-Truth Accuracy (Before vs After Refactoring)",
        fontsize=13,
        fontweight="bold",
        y=0.985,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


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


def plot_fig4_pss_profiles(bench_root: Path, out_path: Path) -> None:
    plt.rcParams.update({"font.size": 9.3, "axes.edgecolor": "#CBD5E1"})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 5.8))

    t_o, r_o, p_o = _load_pss_trace(bench_root / "g1_encff400bqq_t8_abba3" / "baseline.rep2.pss.tsv")
    t_n, r_n, p_n = _load_pss_trace(bench_root / "g1_encff400bqq_t8_abba3" / "candidate.rep2.pss.tsv")
    t_w, r_w, p_w = _load_pss_trace(bench_root / "g1_encff400bqq_t8_warm_abba3" / "candidate.rep2.pss.tsv")

    ax1.plot(t_o, p_o, color=C_ORIG, linewidth=2.0, label=f"Original PSS ({t_o[-1]:.1f}m, peak {p_o.max():.2f} GiB)")
    ax1.plot(t_o, r_o, color=C_ORIG, linewidth=1.2, linestyle="--", alpha=0.65, label=f"Original RSS (peak {r_o.max():.2f} GiB)")
    ax1.plot(t_n, p_n, color=C_NEXT, linewidth=2.2, label=f"Refactored Cold .iqidx PSS ({t_n[-1]:.2f}m, peak {p_n.max():.2f} GiB)")
    ax1.plot(t_w, p_w, color=C_ACCENT, linewidth=2.2, label=f"Refactored Warm .iqidx PSS ({t_w[-1]:.2f}m, peak {p_w.max():.2f} GiB)")
    ax1.set_xlabel("Elapsed Wall-Clock Time (minutes)")
    ax1.set_ylabel("Process-Tree Memory (GiB)")
    ax1.set_ylim(-0.4, 15.2)
    ax1.set_title("A. Tier 1 ENCFF400BQQ (-t 8) Time-Series PSS / RSS Trajectory", fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.legend(loc="upper right", frameon=True, framealpha=0.95, fontsize=8.2)

    t_so, r_so, p_so = _load_pss_trace(bench_root / "g3_shrbmx_rep2_v46_t32" / "baseline.rep1.pss.tsv")
    t_sn, r_sn, p_sn = _load_pss_trace(bench_root / "g3_shrbmx_rep2_v46_t32" / "candidate.rep1.pss.tsv")
    ax2.plot(t_so, r_so, color=C_ORIG, linewidth=1.6, linestyle="--", alpha=0.75, label=f"Original RSS (peak {r_so.max():.2f} GiB)")
    ax2.plot(t_so, p_so, color=C_ORIG, linewidth=2.0, label=f"Original PSS ({t_so[-1]:.1f}m, peak {p_so.max():.2f} GiB)")
    ax2.plot(t_sn, r_sn, color=C_NEXT, linewidth=1.6, linestyle="--", alpha=0.85, label=f"Refactored RSS (peak {r_sn.max():.2f} GiB, -65.1%)")
    ax2.plot(t_sn, p_sn, color=C_ACCENT, linewidth=2.2, label=f"Refactored PSS ({t_sn[-1]:.1f}m, peak {p_sn.max():.2f} GiB, 4.71x faster)")
    ax2.set_xlabel("Elapsed Wall-Clock Time (minutes)")
    ax2.set_ylabel("Process-Tree Memory (GiB)")
    ax2.set_ylim(-1.5, 69.0)
    ax2.set_title("B. Tier 3 shRBMX_rep2 (20.44M ONT + v46, -t 32) Trajectory", fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.legend(loc="upper right", frameon=True, framealpha=0.95, fontsize=8.2)

    fig.suptitle("Process-Tree Memory Trajectories Over Execution Time (/proc/<pid>/smaps_rollup)", fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench-root", type=Path, required=True)
    ap.add_argument("--scatter-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = _load_all(args.bench_root)

    plot_fig1_speedup(rows, args.out_dir / "fig1_speedup_and_walltime.png")
    plot_fig2_memory_and_scaling(rows, args.out_dir / "fig2_memory_and_scaling.png")
    plot_fig3_correlation(rows, args.scatter_dir, args.bench_root, args.out_dir / "fig3_correlation_and_accuracy.png")
    plot_fig4_pss_profiles(args.bench_root, args.out_dir / "fig4_pss_memory_profiles.png")
    print(f"Generated 4 figures in {args.out_dir}")


if __name__ == "__main__":
    main()
