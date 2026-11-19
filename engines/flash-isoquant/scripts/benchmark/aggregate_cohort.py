#!/usr/bin/env python3
"""Aggregate all Slurm benchmark results across all_samples_manifest.json into JSON and TSV tables."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional


def _safe_mean(vals: List[float]) -> Optional[float]:
    return sum(vals) / len(vals) if vals else None


def _safe_min(vals: List[float]) -> Optional[float]:
    return min(vals) if vals else None


def _safe_max(vals: List[float]) -> Optional[float]:
    return max(vals) if vals else None


def _safe_median(vals: List[float]) -> Optional[float]:
    return statistics.median(vals) if vals else None


def _extract_acc_metrics(acc_docs: List[dict]) -> dict:
    gene_rs: List[float] = []
    gene_exacts: List[float] = []
    tx_rs: List[float] = []
    tx_exacts: List[float] = []
    disc_tx_rs: List[float] = []
    disc_exon_rs: List[float] = []
    tm_js: List[float] = []
    ext_js: List[float] = []
    ri_genes: List[float] = []
    ri_fulls: List[float] = []

    for doc in acc_docs:
        gc = doc.get("gene_counts", {})
        if gc.get("present") and gc.get("pearson_r") is not None:
            gene_rs.append(float(gc["pearson_r"]))
            gene_exacts.append(float(gc["exact_match_fraction"]))
        tc = doc.get("transcript_counts", {})
        if tc.get("present") and tc.get("pearson_r") is not None:
            tx_rs.append(float(tc["pearson_r"]))
            tx_exacts.append(float(tc["exact_match_fraction"]))
        dt = doc.get("discovered_transcript_counts", {})
        if dt.get("present") and dt.get("pearson_r") is not None:
            disc_tx_rs.append(float(dt["pearson_r"]))
        de = doc.get("discovered_transcript_by_exons", {})
        if de.get("present") and de.get("pearson_r") is not None:
            disc_exon_rs.append(float(de["pearson_r"]))
        tm = doc.get("transcript_models_gtf", {})
        if tm.get("present") and tm.get("jaccard") is not None:
            tm_js.append(float(tm["jaccard"]))
        ea = doc.get("extended_annotation_gtf", {})
        if ea.get("present") and ea.get("jaccard") is not None:
            ext_js.append(float(ea["jaccard"]))
        ri = doc.get("read_info", {})
        if ri.get("present") and ri.get("gene_agreement") is not None:
            ri_genes.append(float(ri["gene_agreement"]))
            ri_fulls.append(float(ri["full_assignment_agreement"]))

    # Primary expression correlation across annotated or de novo runs
    primary_tx_r = _safe_min(tx_rs) if tx_rs else (_safe_min(disc_exon_rs) if disc_exon_rs else _safe_min(disc_tx_rs))
    primary_gene_r = _safe_min(gene_rs) if gene_rs else _safe_min(disc_tx_rs)

    return {
        "n_acc_files": len(acc_docs),
        "gene_pearson_r": _safe_min(gene_rs),
        "gene_exact_fraction": _safe_min(gene_exacts),
        "transcript_pearson_r": _safe_min(tx_rs),
        "transcript_exact_fraction": _safe_min(tx_exacts),
        "discovered_tx_pearson_r": _safe_min(disc_tx_rs),
        "discovered_exon_pearson_r": _safe_min(disc_exon_rs),
        "primary_gene_or_disc_r": primary_gene_r,
        "primary_tx_or_disc_r": primary_tx_r,
        "transcript_models_jaccard": _safe_min(tm_js),
        "extended_annotation_jaccard": _safe_min(ext_js),
        "read_info_gene_agreement": _safe_min(ri_genes),
        "read_info_full_agreement": _safe_min(ri_fulls),
    }


def _extract_rusage(runs: List[dict], label: str) -> tuple[Optional[float], Optional[float]]:
    sub = [r for r in runs if r.get("label") == label and r.get("exit_code") == 0]
    if not sub:
        return None, None
    u = statistics.median(float(r["user_seconds"]) for r in sub if r.get("user_seconds") is not None)
    s = statistics.median(float(r["sys_seconds"]) for r in sub if r.get("sys_seconds") is not None)
    return round(u, 4), round(s, 4)


def aggregate(manifest_path: Path, bench_root: Path, strict: bool = False) -> dict:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows: List[dict] = []
    missing_slugs: List[str] = []
    failed_slugs: List[str] = []

    for task in manifest:
        slug = task["slug"]
        tdir = bench_root / slug
        bpath = tdir / "benchmark_summary.json"
        if not bpath.is_file():
            missing_slugs.append(slug)
            continue

        bdoc = json.loads(bpath.read_text(encoding="utf-8"))
        summary = bdoc.get("summary", {})
        runs = bdoc.get("runs", [])
        base = summary.get("baseline")
        cand = summary.get("candidate")
        if not cand or cand.get("failures", 0) > 0 or (base and base.get("failures", 0) > 0):
            failed_slugs.append(slug)
            continue

        acc_docs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(tdir.glob("accuracy*.json"))]
        acc = _extract_acc_metrics(acc_docs)

        eval_path = tdir / f"eval_{task['sample_id']}.json"
        sirv_eval = json.loads(eval_path.read_text(encoding="utf-8")) if eval_path.is_file() else None
        cage_path = tdir / "cage_summary.json"
        cage_info = json.loads(cage_path.read_text(encoding="utf-8")) if cage_path.is_file() else None
        vis_path = tdir / "visualize_summary.json"
        vis_info = json.loads(vis_path.read_text(encoding="utf-8")) if vis_path.is_file() else None

        next_wall = float(cand["wall_median_s"])
        next_pss_kib = float(cand["pss_median_kib"])
        next_rss_kib = float(cand["rss_median_kib"])
        next_user, next_sys = _extract_rusage(runs, "candidate")

        is_paired_engine_cmp = base is not None and task["mode"] != "resume_vis"
        orig_wall = float(base["wall_median_s"]) if base else None
        orig_pss_kib = float(base["pss_median_kib"]) if base else None
        orig_rss_kib = float(base["rss_median_kib"]) if base else None
        orig_user, orig_sys = _extract_rusage(runs, "baseline") if base else (None, None)

        speedup = (orig_wall / next_wall) if (orig_wall is not None and next_wall > 0) else None
        pss_savings = ((orig_pss_kib - next_pss_kib) / orig_pss_kib * 100.0) if (orig_pss_kib and orig_pss_kib > 0) else None
        rss_savings = ((orig_rss_kib - next_rss_kib) / orig_rss_kib * 100.0) if (orig_rss_kib and orig_rss_kib > 0) else None

        rows.append(
            {
                "slug": slug,
                "gate": task["gate"],
                "protocol": task["protocol"],
                "species": task["species"],
                "cell_line": task["cell_line"],
                "sample_id": task["sample_id"],
                "mode": task["mode"],
                "threads": task["threads"],
                "repeats": task["repeats"],
                "size_gib": task["size_gib"],
                "is_paired_engine_cmp": is_paired_engine_cmp,
                "orig_wall_s": round(orig_wall, 4) if orig_wall is not None else None,
                "next_wall_s": round(next_wall, 4),
                "speedup_x": round(speedup, 4) if speedup is not None else None,
                "orig_user_s": orig_user,
                "next_user_s": next_user,
                "orig_sys_s": orig_sys,
                "next_sys_s": next_sys,
                "orig_pss_gib": round(orig_pss_kib / 1048576.0, 4) if orig_pss_kib is not None else None,
                "next_pss_gib": round(next_pss_kib / 1048576.0, 4),
                "pss_delta_pct": round(pss_savings, 2) if pss_savings is not None else None,
                "orig_rss_gib": round(orig_rss_kib / 1048576.0, 4) if orig_rss_kib is not None else None,
                "next_rss_gib": round(next_rss_kib / 1048576.0, 4),
                "rss_delta_pct": round(rss_savings, 2) if rss_savings is not None else None,
                **acc,
                "sirv_eval": sirv_eval,
                "cage_info": cage_info,
                "vis_info": vis_info,
                "notes": task.get("notes", ""),
            }
        )

    if strict and (missing_slugs or failed_slugs):
        raise RuntimeError(f"Strict check failed: missing={missing_slugs}, failed={failed_slugs}")

    def _summarize_group(key_fn) -> List[dict]:
        grouped: Dict[str, List[dict]] = defaultdict(list)
        for r in rows:
            grouped[key_fn(r)].append(r)
        out: List[dict] = []
        for gname, items in grouped.items():
            paired = [x for x in items if x["is_paired_engine_cmp"]]
            spds = [x["speedup_x"] for x in paired if x["speedup_x"] is not None]
            pss_d = [x["pss_delta_pct"] for x in paired if x["pss_delta_pct"] is not None]
            rss_d = [x["rss_delta_pct"] for x in paired if x["rss_delta_pct"] is not None]
            ow = [x["orig_wall_s"] for x in paired if x["orig_wall_s"] is not None]
            nw = [x["next_wall_s"] for x in paired]
            op = [x["orig_pss_gib"] for x in paired if x["orig_pss_gib"] is not None]
            np_ = [x["next_pss_gib"] for x in paired]
            or_ = [x["orig_rss_gib"] for x in paired if x["orig_rss_gib"] is not None]
            nr_ = [x["next_rss_gib"] for x in paired]
            rs = [x["primary_tx_or_disc_r"] for x in paired if x["primary_tx_or_disc_r"] is not None]
            grs = [x["primary_gene_or_disc_r"] for x in paired if x["primary_gene_or_disc_r"] is not None]
            ris = [x["read_info_full_agreement"] for x in paired if x["read_info_full_agreement"] is not None]
            tmj = [x["transcript_models_jaccard"] for x in paired if x["transcript_models_jaccard"] is not None]
            exj = [x["extended_annotation_jaccard"] for x in paired if x["extended_annotation_jaccard"] is not None]
            out.append(
                {
                    "group": gname,
                    "n_total": len(items),
                    "n_paired": len(paired),
                    "orig_wall_median_s": round(_safe_median(ow), 2) if ow else None,
                    "next_wall_median_s": round(_safe_median(nw), 2) if nw else None,
                    "orig_wall_sum_s": round(sum(ow), 1) if ow else None,
                    "next_wall_sum_s": round(sum(nw), 1) if nw else None,
                    "speedup_median_x": round(_safe_median(spds), 2) if spds else None,
                    "speedup_min_x": round(_safe_min(spds), 2) if spds else None,
                    "speedup_max_x": round(_safe_max(spds), 2) if spds else None,
                    "orig_pss_median_gib": round(_safe_median(op), 2) if op else None,
                    "next_pss_median_gib": round(_safe_median(np_), 2) if np_ else None,
                    "pss_delta_median_pct": round(_safe_median(pss_d), 1) if pss_d else None,
                    "orig_rss_median_gib": round(_safe_median(or_), 2) if or_ else None,
                    "next_rss_median_gib": round(_safe_median(nr_), 2) if nr_ else None,
                    "rss_delta_median_pct": round(_safe_median(rss_d), 1) if rss_d else None,
                    "gene_r_min": round(_safe_min(grs), 8) if grs else None,
                    "tx_r_min": round(_safe_min(rs), 8) if rs else None,
                    "read_info_min": round(_safe_min(ris), 6) if ris else None,
                    "tx_models_jaccard_min": round(_safe_min(tmj), 6) if tmj else None,
                    "ext_gtf_jaccard_min": round(_safe_min(exj), 6) if exj else None,
                }
            )
        return out

    paired_all = [r for r in rows if r["is_paired_engine_cmp"]]
    global_summary = {
        "manifest_total": len(manifest),
        "completed_tasks": len(rows),
        "paired_engine_tasks": len(paired_all),
        "missing_count": len(missing_slugs),
        "failed_count": len(failed_slugs),
        "missing_slugs": missing_slugs,
        "failed_slugs": failed_slugs,
        "total_orig_wall_hours": round(sum(r["orig_wall_s"] for r in paired_all if r["orig_wall_s"]) / 3600.0, 2),
        "total_next_wall_hours": round(sum(r["next_wall_s"] for r in paired_all) / 3600.0, 2),
        "overall_speedup_median_x": round(_safe_median([r["speedup_x"] for r in paired_all if r["speedup_x"]]), 2) if paired_all else None,
        "overall_speedup_min_x": round(_safe_min([r["speedup_x"] for r in paired_all if r["speedup_x"]]), 2) if paired_all else None,
        "overall_speedup_max_x": round(_safe_max([r["speedup_x"] for r in paired_all if r["speedup_x"]]), 2) if paired_all else None,
        "overall_pss_delta_median_pct": round(_safe_median([r["pss_delta_pct"] for r in paired_all if r["pss_delta_pct"] is not None]), 1) if paired_all else None,
        "overall_rss_delta_median_pct": round(_safe_median([r["rss_delta_pct"] for r in paired_all if r["rss_delta_pct"] is not None]), 1) if paired_all else None,
        "overall_tx_r_min": round(_safe_min([r["primary_tx_or_disc_r"] for r in paired_all if r["primary_tx_or_disc_r"] is not None]), 8) if paired_all else None,
        "overall_gene_r_min": round(_safe_min([r["primary_gene_or_disc_r"] for r in paired_all if r["primary_gene_or_disc_r"] is not None]), 8) if paired_all else None,
        "overall_read_info_min": round(_safe_min([r["read_info_full_agreement"] for r in paired_all if r["read_info_full_agreement"] is not None]), 6) if paired_all else None,
    }

    return {
        "global": global_summary,
        "by_gate": _summarize_group(lambda r: r["gate"]),
        "by_protocol": _summarize_group(lambda r: r["protocol"]),
        "by_species": _summarize_group(lambda r: r["species"]),
        "rows": rows,
    }


def _write_tsv(path: Path, records: List[dict], cols: List[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in records:
            w.writerow([r.get(c) for c in cols])


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate full-cohort Slurm benchmark metrics.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--bench-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--strict", action="store_true", default=False)
    args = parser.parse_args()

    out_dir = args.out_dir or args.bench_root
    out_dir.mkdir(parents=True, exist_ok=True)
    res = aggregate(args.manifest, args.bench_root, strict=args.strict)

    (out_dir / "cohort_all_results.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")

    row_cols = [
        "slug", "gate", "protocol", "species", "cell_line", "sample_id", "mode", "threads",
        "size_gib", "is_paired_engine_cmp", "orig_wall_s", "next_wall_s", "speedup_x",
        "orig_pss_gib", "next_pss_gib", "pss_delta_pct", "orig_rss_gib", "next_rss_gib", "rss_delta_pct",
        "primary_gene_or_disc_r", "primary_tx_or_disc_r", "read_info_full_agreement",
        "transcript_models_jaccard", "extended_annotation_jaccard",
    ]
    _write_tsv(out_dir / "cohort_all_results.tsv", res["rows"], row_cols)

    grp_cols = [
        "group", "n_total", "n_paired", "orig_wall_median_s", "next_wall_median_s",
        "orig_wall_sum_s", "next_wall_sum_s", "speedup_median_x", "speedup_min_x", "speedup_max_x",
        "orig_pss_median_gib", "next_pss_median_gib", "pss_delta_median_pct",
        "orig_rss_median_gib", "next_rss_median_gib", "rss_delta_median_pct",
        "gene_r_min", "tx_r_min", "read_info_min", "tx_models_jaccard_min", "ext_gtf_jaccard_min",
    ]
    _write_tsv(out_dir / "cohort_by_protocol.tsv", res["by_protocol"], grp_cols)
    _write_tsv(out_dir / "cohort_by_gate.tsv", res["by_gate"], grp_cols)

    print(json.dumps(res["global"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
