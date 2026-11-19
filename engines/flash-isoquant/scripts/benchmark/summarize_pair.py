#!/usr/bin/env python3
"""Compare original vs isoquant-next output directories and emit quantitative accuracy metrics."""

from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

from compare_outputs import (
    _extract_gtf_attr,
    _load_gtf_structures,
    compare_count_tables,
    compare_gtf_models,
    compare_read_assignments,
)


def _load_single_col_counts(path: Path) -> Dict[str, float]:
    out: Dict[str, float] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#") or line.startswith("__"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                out[parts[0]] = float(parts[1])
            except ValueError:
                continue
    return out


def _pearson(xs: List[float], ys: List[float]) -> float:
    n = len(xs)
    if n < 2:
        return 1.0 if xs == ys else 0.0
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = sum((x - mx) ** 2 for x in xs)
    dy = sum((y - my) ** 2 for y in ys)
    if dx == 0.0 and dy == 0.0:
        return 1.0
    if dx == 0.0 or dy == 0.0:
        return 0.0
    return num / math.sqrt(dx * dy)


def _quant_metrics(ref_path: Path, cand_path: Path) -> dict:
    if not ref_path.exists() or not cand_path.exists():
        return {"present": False, "ref_exists": ref_path.exists(), "cand_exists": cand_path.exists()}
    comp = compare_count_tables(ref_path, cand_path)
    ref_m = _load_single_col_counts(ref_path)
    cand_m = _load_single_col_counts(cand_path)
    keys = sorted(set(ref_m) | set(cand_m))
    xs = [ref_m.get(k, 0.0) for k in keys]
    ys = [cand_m.get(k, 0.0) for k in keys]
    exact = sum(1 for a, b in zip(xs, ys) if abs(a - b) <= 1e-6)
    max_diff = max((abs(a - b) for a, b in zip(xs, ys)), default=0.0)
    return {
        "present": True,
        "compare_passed": comp.passed,
        "compare_details": comp.details,
        "num_features": len(keys),
        "exact_match_fraction": round(exact / max(1, len(keys)), 6),
        "max_abs_diff": round(max_diff, 6),
        "pearson_r": round(_pearson(xs, ys), 8),
        "ref_sum": round(sum(xs), 4),
        "cand_sum": round(sum(ys), 4),
    }


def _load_tx_signatures(gtf_path: Path) -> Dict[str, Tuple[str, str, Tuple[Tuple[int, int], ...]]]:
    tx_exons: Dict[str, Tuple[str, str, List[Tuple[int, int]]]] = {}
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    with opener(gtf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            tid = _extract_gtf_attr(cols[8], "transcript_id")
            if not tid:
                continue
            exon = (int(cols[3]), int(cols[4]))
            entry = tx_exons.get(tid)
            if entry is None:
                tx_exons[tid] = (cols[0], cols[6], [exon])
            else:
                entry[2].append(exon)
    return {
        tid: (chrom, strand, tuple(sorted(exons)))
        for tid, (chrom, strand, exons) in tx_exons.items()
    }


def _discovered_exon_count_metrics(
    ref_gtf: Path, ref_counts: Path, cand_gtf: Path, cand_counts: Path
) -> dict:
    if not (ref_gtf.exists() and ref_counts.exists() and cand_gtf.exists() and cand_counts.exists()):
        return {"present": False}
    ref_sigs = _load_tx_signatures(ref_gtf)
    cand_sigs = _load_tx_signatures(cand_gtf)
    ref_c = _load_single_col_counts(ref_counts)
    cand_c = _load_single_col_counts(cand_counts)

    ref_by_sig: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], float] = {}
    for tid, val in ref_c.items():
        sig = ref_sigs.get(tid)
        if sig is not None:
            ref_by_sig[sig] = ref_by_sig.get(sig, 0.0) + val

    cand_by_sig: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], float] = {}
    for tid, val in cand_c.items():
        sig = cand_sigs.get(tid)
        if sig is not None:
            cand_by_sig[sig] = cand_by_sig.get(sig, 0.0) + val

    keys = list(set(ref_by_sig) | set(cand_by_sig))
    xs = [ref_by_sig.get(k, 0.0) for k in keys]
    ys = [cand_by_sig.get(k, 0.0) for k in keys]
    exact = sum(1 for a, b in zip(xs, ys) if abs(a - b) <= 1e-6)
    return {
        "present": True,
        "num_chains": len(keys),
        "exact_match_fraction": round(exact / max(1, len(keys)), 6),
        "pearson_r": round(_pearson(xs, ys), 8),
        "ref_sum": round(sum(xs), 4),
        "cand_sum": round(sum(ys), 4),
    }


def _gtf_metrics(ref_path: Path, cand_path: Path) -> dict:
    if not ref_path.exists() or not cand_path.exists():
        return {"present": False, "ref_exists": ref_path.exists(), "cand_exists": cand_path.exists()}
    comp = compare_gtf_models(ref_path, cand_path)
    ref_s = set(_load_gtf_structures(ref_path))
    cand_s = set(_load_gtf_structures(cand_path))
    inter = len(ref_s & cand_s)
    union = len(ref_s | cand_s)
    return {
        "present": True,
        "compare_passed": comp.passed,
        "compare_details": comp.details,
        "ref_models": len(ref_s),
        "cand_models": len(cand_s),
        "shared_models": inter,
        "jaccard": round(inter / max(1, union), 6),
    }


def _assign_metrics(ref_path: Path, cand_path: Path) -> dict:
    if not ref_path.exists() or not cand_path.exists():
        return {"present": False, "ref_exists": ref_path.exists(), "cand_exists": cand_path.exists()}
    # Stream sorted-key map cleanly without huge memory blowup if small enough; compare exact tuple (read_id, gene_id, tx_id, type)
    def _load_keyed(path: Path) -> Dict[Tuple[str, str], Tuple[str, str, str]]:
        opener = gzip.open if path.suffix == ".gz" else open
        out: Dict[Tuple[str, str], Tuple[str, str, str]] = {}
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 7:
                    key = (cols[0], cols[1])
                    out[key] = (cols[3], cols[4], cols[5])
        return out

    ref_map = _load_keyed(ref_path)
    cand_map = _load_keyed(cand_path)
    all_keys = set(ref_map) | set(cand_map)
    exact_tuple = sum(1 for k in all_keys if ref_map.get(k) == cand_map.get(k))
    same_gene = sum(1 for k in all_keys if ref_map.get(k, ("",))[0] == cand_map.get(k, ("?",))[0])
    # Run exact compare_read_assignments on small/medium (<2M records) or when equal counts
    exact_full = None
    if len(all_keys) <= 2_500_000:
        exact_full = compare_read_assignments(ref_path, cand_path).passed
    return {
        "present": True,
        "ref_records": len(ref_map),
        "cand_records": len(cand_map),
        "gene_agreement": round(same_gene / max(1, len(all_keys)), 6),
        "full_assignment_agreement": round(exact_tuple / max(1, len(all_keys)), 6),
        "exact_sorted_line_match": exact_full,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize accuracy between reference and candidate IsoQuant runs.")
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--cand-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()

    s_ref = args.ref_dir / args.prefix
    s_cand = args.cand_dir / args.prefix

    payload = {
        "prefix": args.prefix,
        "gene_counts": _quant_metrics(s_ref / f"{args.prefix}.gene_counts.tsv", s_cand / f"{args.prefix}.gene_counts.tsv"),
        "transcript_counts": _quant_metrics(
            s_ref / f"{args.prefix}.transcript_counts.tsv", s_cand / f"{args.prefix}.transcript_counts.tsv"
        ),
        "gene_tpm": _quant_metrics(s_ref / f"{args.prefix}.gene_tpm.tsv", s_cand / f"{args.prefix}.gene_tpm.tsv"),
        "transcript_tpm": _quant_metrics(
            s_ref / f"{args.prefix}.transcript_tpm.tsv", s_cand / f"{args.prefix}.transcript_tpm.tsv"
        ),
        "transcript_models_gtf": _gtf_metrics(
            s_ref / f"{args.prefix}.transcript_models.gtf", s_cand / f"{args.prefix}.transcript_models.gtf"
        ),
        "discovered_transcript_counts": _quant_metrics(
            s_ref / f"{args.prefix}.discovered_transcript_counts.tsv",
            s_cand / f"{args.prefix}.discovered_transcript_counts.tsv",
        ),
        "discovered_transcript_tpm": _quant_metrics(
            s_ref / f"{args.prefix}.discovered_transcript_tpm.tsv",
            s_cand / f"{args.prefix}.discovered_transcript_tpm.tsv",
        ),
        "discovered_gene_counts": _quant_metrics(
            s_ref / f"{args.prefix}.discovered_gene_counts.tsv",
            s_cand / f"{args.prefix}.discovered_gene_counts.tsv",
        ),
        "discovered_transcript_by_exons": _discovered_exon_count_metrics(
            s_ref / f"{args.prefix}.transcript_models.gtf",
            s_ref / f"{args.prefix}.discovered_transcript_counts.tsv",
            s_cand / f"{args.prefix}.transcript_models.gtf",
            s_cand / f"{args.prefix}.discovered_transcript_counts.tsv",
        ),
        "extended_annotation_gtf": _gtf_metrics(
            s_ref / f"{args.prefix}.extended_annotation.gtf", s_cand / f"{args.prefix}.extended_annotation.gtf"
        ),
        "read_info": _assign_metrics(
            s_ref / f"{args.prefix}.read_info.tsv.gz", s_cand / f"{args.prefix}.read_info.tsv.gz"
        ),
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
