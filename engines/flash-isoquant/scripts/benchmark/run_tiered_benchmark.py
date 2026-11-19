#!/usr/bin/env python3
"""Run alternating A-B-B-A benchmarks with process-tree PSS/RSS and wall/user/sys accounting."""

from __future__ import annotations

import argparse
import csv
import json
import os
import statistics
import subprocess
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional


MONITOR_SH = Path(__file__).resolve().with_name("monitor_process_tree.sh")


@dataclass
class RunMetrics:
    label: str
    rep: int
    exit_code: int
    wall_seconds: float
    user_seconds: float
    sys_seconds: float
    peak_tree_pss_kib: int
    peak_tree_rss_kib: int
    max_process_count: int
    pss_trace_path: str
    stdout_log: str
    stderr_log: str


def run_monitored_command(
    cmd: List[str],
    label: str,
    rep: int,
    log_dir: Path,
    interval: float = 0.5,
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
) -> RunMetrics:
    log_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{label}.rep{rep}"
    phase_file = log_dir / f"{stem}.phase"
    trace_tsv = log_dir / f"{stem}.pss.tsv"
    stdout_path = log_dir / f"{stem}.stdout.log"
    stderr_path = log_dir / f"{stem}.stderr.log"
    phase_file.write_text("run\n", encoding="utf-8")

    with open(stdout_path, "w", encoding="utf-8") as out_fh, open(
        stderr_path, "w", encoding="utf-8"
    ) as err_fh:
        t0 = time.monotonic()
        proc = subprocess.Popen(
            cmd,
            stdout=out_fh,
            stderr=err_fh,
            cwd=str(cwd) if cwd else None,
            env=env,
        )
        mon = subprocess.Popen(
            ["bash", str(MONITOR_SH), str(proc.pid), str(phase_file), str(trace_tsv), str(interval)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        _, status, rusage = os.wait4(proc.pid, 0)
        wall = time.monotonic() - t0
        mon.wait(timeout=10)

    exit_code = os.waitstatus_to_exitcode(status)
    peak_pss = 0
    peak_rss = 0
    max_procs = 0
    if trace_tsv.exists():
        with open(trace_tsv, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                peak_pss = max(peak_pss, int(row["tree_pss_kib"]))
                peak_rss = max(peak_rss, int(row["tree_rss_kib"]))
                max_procs = max(max_procs, int(row["process_count"]))

    # Fallback when very short commands finish before the first /proc sample
    if peak_rss == 0 and rusage.ru_maxrss > 0:
        peak_rss = int(rusage.ru_maxrss)
    if peak_pss == 0 and rusage.ru_maxrss > 0:
        peak_pss = int(rusage.ru_maxrss)

    return RunMetrics(
        label=label,
        rep=rep,
        exit_code=exit_code,
        wall_seconds=round(wall, 4),
        user_seconds=round(float(rusage.ru_utime), 4),
        sys_seconds=round(float(rusage.ru_stime), 4),
        peak_tree_pss_kib=peak_pss,
        peak_tree_rss_kib=peak_rss,
        max_process_count=max_procs,
        pss_trace_path=str(trace_tsv),
        stdout_log=str(stdout_path),
        stderr_log=str(stderr_path),
    )


def summarize_by_label(runs: List[RunMetrics]) -> Dict[str, Dict[str, float]]:
    grouped: Dict[str, List[RunMetrics]] = {}
    for r in runs:
        grouped.setdefault(r.label, []).append(r)

    summary: Dict[str, Dict[str, float]] = {}
    for label, items in grouped.items():
        walls = [x.wall_seconds for x in items]
        pss_vals = [x.peak_tree_pss_kib for x in items]
        rss_vals = [x.peak_tree_rss_kib for x in items]
        summary[label] = {
            "count": float(len(items)),
            "failures": float(sum(1 for x in items if x.exit_code != 0)),
            "wall_median_s": statistics.median(walls),
            "wall_min_s": min(walls),
            "wall_max_s": max(walls),
            "pss_median_kib": statistics.median(pss_vals),
            "pss_min_kib": float(min(pss_vals)),
            "pss_max_kib": float(max(pss_vals)),
            "rss_median_kib": statistics.median(rss_vals),
        }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Run alternating monitored benchmark commands.")
    parser.add_argument("--baseline-cmd", help="Shell command string for baseline A")
    parser.add_argument("--candidate-cmd", help="Shell command string for candidate B")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--interval", type=float, default=0.5)
    parser.add_argument("--log-dir", type=Path, required=True)
    args = parser.parse_args()

    if not args.baseline_cmd and not args.candidate_cmd:
        parser.error("Provide --baseline-cmd and/or --candidate-cmd.")

    runs: List[RunMetrics] = []
    for rep in range(1, args.repeats + 1):
        order: List[tuple[str, str]] = []
        if args.baseline_cmd and args.candidate_cmd:
            pair = [("baseline", args.baseline_cmd), ("candidate", args.candidate_cmd)]
            order = pair if rep % 2 == 1 else list(reversed(pair))
        elif args.baseline_cmd:
            order = [("baseline", args.baseline_cmd)]
        elif args.candidate_cmd:
            order = [("candidate", args.candidate_cmd)]

        for label, cmd_str in order:
            runs.append(
                run_monitored_command(
                    ["bash", "-lc", cmd_str],
                    label=label,
                    rep=rep,
                    log_dir=args.log_dir,
                    interval=args.interval,
                )
            )

    report = {
        "summary": summarize_by_label(runs),
        "runs": [asdict(r) for r in runs],
    }
    out_json = args.log_dir / "benchmark_summary.json"
    out_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
