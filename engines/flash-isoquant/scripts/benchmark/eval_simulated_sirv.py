#!/usr/bin/env python3
"""Evaluate simulated NanoSim read-assignment and transcript-discovery Recall/Precision."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Dict, Set


def _base_id(tx_id: str) -> str:
    return tx_id.split(".", 1)[0]


def _extract_gtf_tx_info(gtf_path: Path) -> tuple[Set[str], Dict[str, str]]:
    """Return discovered transcript IDs and transcript -> similar_reference_id (if novel)."""
    tx_ids: Set[str] = set()
    sim_refs: Dict[str, str] = {}
    with open(gtf_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            attrs = cols[8]
            tid = ""
            if 'transcript_id "' in attrs:
                tid = attrs.split('transcript_id "', 1)[1].split('"', 1)[0]
            if not tid:
                continue
            tx_ids.add(tid)
            if 'similar_reference_id "' in attrs:
                sim_refs[tid] = attrs.split('similar_reference_id "', 1)[1].split('"', 1)[0]
    return tx_ids, sim_refs


def evaluate_run(run_dir: Path, prefix: str) -> dict:
    sdir = run_dir / prefix
    ri_path = sdir / f"{prefix}.read_info.tsv.gz"
    gtf_path = sdir / f"{prefix}.transcript_models.gtf"
    tc_path = sdir / f"{prefix}.transcript_counts.tsv"

    total_reads = 0
    assigned_unique = 0
    correct_unique = 0
    ambiguous = 0
    inconsistent = 0
    noninformative = 0
    true_tx_read_counts: Dict[str, int] = {}

    with gzip.open(ri_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            read_id = cols[0]
            true_base = read_id.split("_", 1)[0]
            true_tx_read_counts[true_base] = true_tx_read_counts.get(true_base, 0) + 1
            total_reads += 1

            pred_tx = _base_id(cols[5])
            atype = cols[6]
            if atype in ("unique", "unique_minor_difference"):
                assigned_unique += 1
                if pred_tx == true_base:
                    correct_unique += 1
            elif atype == "ambiguous":
                ambiguous += 1
            elif atype == "noninformative":
                noninformative += 1
            else:
                inconsistent += 1

    # Transcript discovery & expressed known-isoform recall/precision
    true_expressed_ge2 = {tx for tx, cnt in true_tx_read_counts.items() if cnt >= 2}
    disc_tx_ids, _ = _extract_gtf_tx_info(gtf_path) if gtf_path.exists() else (set(), {})
    disc_bases = {_base_id(t) for t in disc_tx_ids if not t.startswith("transcript")}

    recovered_ge2 = len(disc_bases & true_expressed_ge2)
    tx_recall_ge2 = recovered_ge2 / max(1, len(true_expressed_ge2))
    known_disc_precision = len(disc_bases & set(true_tx_read_counts)) / max(1, len(disc_tx_ids))

    # Count-table expressed transcripts (count >= 1)
    quant_positives: Set[str] = set()
    if tc_path.exists():
        with open(tc_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if not line or line.startswith("#") or line.startswith("__"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    try:
                        if float(parts[1]) >= 1.0:
                            quant_positives.add(_base_id(parts[0]))
                    except ValueError:
                        pass

    quant_tp = len(quant_positives & set(true_tx_read_counts))
    quant_recall_ge2 = len(quant_positives & true_expressed_ge2) / max(1, len(true_expressed_ge2))
    quant_precision = quant_tp / max(1, len(quant_positives))

    return {
        "total_read_records": total_reads,
        "assigned_unique": assigned_unique,
        "correct_unique_transcript": correct_unique,
        "read_precision": round(correct_unique / max(1, assigned_unique), 6),
        "read_recall": round(correct_unique / max(1, total_reads), 6),
        "ambiguous": ambiguous,
        "inconsistent": inconsistent,
        "noninformative": noninformative,
        "true_expressed_tx_ge2_reads": len(true_expressed_ge2),
        "discovered_total_models": len(disc_tx_ids),
        "discovered_known_models": len(disc_bases),
        "model_recall_ge2_reads": round(tx_recall_ge2, 6),
        "model_known_over_all_precision": round(known_disc_precision, 6),
        "quant_expressed_tx": len(quant_positives),
        "quant_recall_ge2_reads": round(quant_recall_ge2, 6),
        "quant_precision": round(quant_precision, 6),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate NanoSim ground truth accuracy for IsoQuant runs.")
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--cand-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()

    out = {
        "baseline": evaluate_run(args.ref_dir, args.prefix),
        "candidate": evaluate_run(args.cand_dir, args.prefix),
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
