#!/usr/bin/env python3
"""Strict equivalence diff checker between IsoQuant baseline and isoquant_next outputs."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass
class ComparisonResult:
    name: str
    passed: bool
    details: str


def _open_maybe_gz(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def compare_count_tables(
    ref_path: Path, cand_path: Path, float_atol: float = 1e-6
) -> ComparisonResult:
    """Compare TSV count matrices (exact row/col IDs and numerical equality)."""
    if not ref_path.exists() or not cand_path.exists():
        return ComparisonResult(
            name=ref_path.name,
            passed=False,
            details=f"missing file: ref={ref_path.exists()} cand={cand_path.exists()}",
        )

    def load_table(p: Path) -> Tuple[List[str], Dict[str, List[float]]]:
        header: List[str] = []
        rows: Dict[str, List[float]] = {}
        with _open_maybe_gz(p) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if not header:
                    header = parts
                    continue
                rows[parts[0]] = [float(x) for x in parts[1:]]
        return header, rows

    ref_hdr, ref_rows = load_table(ref_path)
    cand_hdr, cand_rows = load_table(cand_path)

    if ref_hdr != cand_hdr:
        return ComparisonResult(
            name=ref_path.name,
            passed=False,
            details=f"header mismatch: ref={ref_hdr} cand={cand_hdr}",
        )
    ref_keys = set(ref_rows)
    cand_keys = set(cand_rows)
    if ref_keys != cand_keys:
        missing = sorted(ref_keys - cand_keys)[:5]
        extra = sorted(cand_keys - ref_keys)[:5]
        return ComparisonResult(
            name=ref_path.name,
            passed=False,
            details=f"key set mismatch: missing={missing}, extra={extra}",
        )

    max_diff = 0.0
    mismatch_key: Optional[str] = None
    for k, rv in ref_rows.items():
        cv = cand_rows[k]
        for a, b in zip(rv, cv):
            d = abs(a - b)
            if math.isnan(d) or d > float_atol:
                if d > max_diff or mismatch_key is None:
                    max_diff = d
                    mismatch_key = k
    if mismatch_key is not None:
        return ComparisonResult(
            name=ref_path.name,
            passed=False,
            details=f"value mismatch at {mismatch_key}: max_abs_diff={max_diff:.3e} > {float_atol}",
        )
    return ComparisonResult(
        name=ref_path.name,
        passed=True,
        details=f"exact/tol match across {len(ref_rows)} features",
    )


def _extract_gtf_attr(attrs: str, key: str) -> Optional[str]:
    marker = f'{key} "'
    idx = attrs.find(marker)
    if idx < 0:
        return None
    start = idx + len(marker)
    end = attrs.find('"', start)
    if end < 0:
        return None
    return attrs[start:end]


def _load_gtf_structures(
    path: Path,
) -> Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], int]:
    tx_exons: Dict[str, Tuple[str, str, List[Tuple[int, int]]]] = {}
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            tx_id = _extract_gtf_attr(cols[8], "transcript_id")
            if tx_id is None:
                continue
            exon = (int(cols[3]), int(cols[4]))
            entry = tx_exons.get(tx_id)
            if entry is None:
                tx_exons[tx_id] = (cols[0], cols[6], [exon])
            else:
                entry[2].append(exon)

    sig_counts: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], int] = {}
    for chrom, strand, exons in tx_exons.values():
        exons_sorted = tuple(sorted(exons))
        sig = (chrom, strand, exons_sorted)
        sig_counts[sig] = sig_counts.get(sig, 0) + 1
    return sig_counts


def compare_gtf_models(
    ref_path: Path, cand_path: Path, coord_tol: int = 0
) -> ComparisonResult:
    """Compare transcript exon-chain structures across GTF files."""
    if not ref_path.exists() or not cand_path.exists():
        return ComparisonResult(
            name=ref_path.name,
            passed=False,
            details=f"missing file: ref={ref_path.exists()} cand={cand_path.exists()}",
        )
    ref_sigs = _load_gtf_structures(ref_path)
    cand_sigs = _load_gtf_structures(cand_path)
    if coord_tol == 0:
        if ref_sigs == cand_sigs:
            return ComparisonResult(
                name=ref_path.name,
                passed=True,
                details=f"exact exon-chain match ({len(ref_sigs)} transcripts)",
            )
        ref_set = set(ref_sigs)
        cand_set = set(cand_sigs)
        return ComparisonResult(
            name=ref_path.name,
            passed=False,
            details=(
                f"transcript structure mismatch: ref={len(ref_sigs)} cand={len(cand_sigs)} "
                f"ref_only={len(ref_set - cand_set)} cand_only={len(cand_set - ref_set)}"
            ),
        )
    # Approximate boundary matching when coord_tol > 0
    if len(ref_sigs) != len(cand_sigs):
        return ComparisonResult(
            name=ref_path.name,
            passed=False,
            details=f"transcript count mismatch: ref={len(ref_sigs)} cand={len(cand_sigs)}",
        )
    ref_sorted = sorted(ref_sigs.keys())
    cand_sorted = sorted(cand_sigs.keys())
    for (rc, rs, re), (cc, cs, ce) in zip(ref_sorted, cand_sorted):
        if rc != cc or rs != cs or len(re) != len(ce):
            return ComparisonResult(
                name=ref_path.name,
                passed=False,
                details=f"intron-chain mismatch between {(rc, rs, len(re))} and {(cc, cs, len(ce))}",
            )
        for idx, ((ra, rb), (ca, cb)) in enumerate(zip(re, ce)):
            if idx == 0 and abs(ra - ca) > coord_tol:
                return ComparisonResult(
                    name=ref_path.name,
                    passed=False,
                    details=f"TSS boundary shift {abs(ra - ca)} > {coord_tol}",
                )
            if idx == len(re) - 1 and abs(rb - cb) > coord_tol:
                return ComparisonResult(
                    name=ref_path.name,
                    passed=False,
                    details=f"TES boundary shift {abs(rb - cb)} > {coord_tol}",
                )
            if idx > 0 and ra != ca:
                return ComparisonResult(
                    name=ref_path.name,
                    passed=False,
                    details=f"internal splice acceptor mismatch {ra} != {ca}",
                )
            if idx < len(re) - 1 and rb != cb:
                return ComparisonResult(
                    name=ref_path.name,
                    passed=False,
                    details=f"internal splice donor mismatch {rb} != {cb}",
                )
    return ComparisonResult(
        name=ref_path.name,
        passed=True,
        details=f"intron-chain match within coord_tol={coord_tol} ({len(ref_sigs)} transcripts)",
    )


def compare_read_assignments(ref_path: Path, cand_path: Path) -> ComparisonResult:
    """Compare read assignment TSV(.gz) rows keyed by (read_id, chr, gene_id, isoform_id, assignment_type)."""
    if not ref_path.exists() or not cand_path.exists():
        return ComparisonResult(
            name=ref_path.name,
            passed=False,
            details=f"missing file: ref={ref_path.exists()} cand={cand_path.exists()}",
        )

    def load_rows(p: Path) -> List[str]:
        lines: List[str] = []
        with _open_maybe_gz(p) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                lines.append(line.rstrip("\n"))
        return sorted(lines)

    ref_lines = load_rows(ref_path)
    cand_lines = load_rows(cand_path)
    if ref_lines == cand_lines:
        return ComparisonResult(
            name=ref_path.name,
            passed=True,
            details=f"exact sorted line match ({len(ref_lines)} records)",
        )
    return ComparisonResult(
        name=ref_path.name,
        passed=False,
        details=f"assignment mismatch: ref_lines={len(ref_lines)} cand_lines={len(cand_lines)}",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare IsoQuant baseline and refactored outputs.")
    parser.add_argument("--ref-counts", type=Path)
    parser.add_argument("--cand-counts", type=Path)
    parser.add_argument("--ref-gtf", type=Path)
    parser.add_argument("--cand-gtf", type=Path)
    parser.add_argument("--ref-assign", type=Path)
    parser.add_argument("--cand-assign", type=Path)
    parser.add_argument("--float-atol", type=float, default=1e-6)
    parser.add_argument("--coord-tol", type=int, default=0)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    results: List[ComparisonResult] = []
    if args.ref_counts and args.cand_counts:
        results.append(compare_count_tables(args.ref_counts, args.cand_counts, args.float_atol))
    if args.ref_gtf and args.cand_gtf:
        results.append(compare_gtf_models(args.ref_gtf, args.cand_gtf, args.coord_tol))
    if args.ref_assign and args.cand_assign:
        results.append(compare_read_assignments(args.ref_assign, args.cand_assign))

    if not results:
        parser.error("Provide at least one ref/cand pair to compare.")

    payload = [asdict(r) for r in results]
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    else:
        print(json.dumps(payload, indent=2))

    return 0 if all(r.passed for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
