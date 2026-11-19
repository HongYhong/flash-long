#!/usr/bin/env python3
"""A-B-B-A paired benchmark runner for native vs bambu_next (strict no-/tmp & PSS monitor)."""
import argparse
import json
import os
import statistics
import subprocess
import time
from pathlib import Path

WORK_ROOT = Path("/gpfs/group/data_sci_group/weiwu_group/hongyanhong/bambu_next_work")
MONITOR_SH = Path("/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_next/scripts/benchmark/monitor_process_tree.sh")
CANDIDATE_LIB = Path("/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_next/.r_lib")
CONDA_RSCRIPT = "/gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/envs/bambu/bin/Rscript"


def parse_pss_tsv(tsv_path: Path):
    peak_pss_kib = 0
    peak_rss_kib = 0
    samples = 0
    with open(tsv_path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            rss = int(parts[2])
            pss = int(parts[3])
            peak_rss_kib = max(peak_rss_kib, rss)
            peak_pss_kib = max(peak_pss_kib, pss)
            samples += 1
    return {
        "samples": samples,
        "peak_pss_gib": peak_pss_kib / (1024.0 * 1024.0),
        "peak_rss_gib": peak_rss_kib / (1024.0 * 1024.0),
    }


def run_single_trial(case_id: str, impl: str, rep_idx: int, r_script_path: Path, run_dir: Path, extra_args: list[str]):
    trial_dir = run_dir / f"{case_id}_{impl}_rep{rep_idx}"
    tmp_dir = trial_dir / "tmp"
    out_dir = trial_dir / "out"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["TMPDIR"] = str(tmp_dir)
    env["TMP"] = str(tmp_dir)
    env["TEMP"] = str(tmp_dir)
    env["R_SESSION_TMPDIR"] = str(tmp_dir)
    env["XDG_CACHE_HOME"] = str(trial_dir / "cache")
    if impl == "candidate":
        env["R_LIBS_USER"] = str(CANDIDATE_LIB)
    else:
        env.pop("R_LIBS_USER", None)

    phase_file = trial_dir / "phase.txt"
    phase_file.write_text("init\n")
    pss_tsv = trial_dir / "pss_trace.tsv"
    stdout_log = trial_dir / "stdout.log"
    stderr_log = trial_dir / "stderr.log"
    meta_json = trial_dir / "meta.json"

    cmd = [CONDA_RSCRIPT, str(r_script_path), "--impl", impl, "--out-dir", str(out_dir), "--phase-file", str(phase_file)] + extra_args
    t0 = time.monotonic()
    with open(stdout_log, "w") as fout, open(stderr_log, "w") as ferr:
        proc = subprocess.Popen(cmd, env=env, stdout=fout, stderr=ferr, cwd=str(trial_dir))
        mon = subprocess.Popen(["bash", str(MONITOR_SH), str(proc.pid), str(phase_file), str(pss_tsv), "1"])
        rc = proc.wait()
        mon.wait()
    wall_sec = time.monotonic() - t0

    pss_stats = parse_pss_tsv(pss_tsv) if pss_tsv.exists() else {"samples": 0, "peak_pss_gib": 0.0, "peak_rss_gib": 0.0}
    res = {
        "case_id": case_id,
        "impl": impl,
        "rep": rep_idx,
        "rc": rc,
        "wall_seconds": wall_sec,
        **pss_stats,
        "out_dir": str(out_dir),
    }
    meta_json.write_text(json.dumps(res, indent=2))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-id", required=True)
    ap.add_argument("--r-script", required=True)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("extra_args", nargs=argparse.REMAINDER)
    args = ap.parse_args()

    extra = [x for x in args.extra_args if x != "--"]
    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    r_script = Path(args.r_script)

    schedule = []
    for r in range(1, args.rounds + 1):
        # Strict A-B-B-A / alternating order per pair
        order = ("native", "candidate") if (r % 2 == 1) else ("candidate", "native")
        for impl in order:
            schedule.append((r, impl))

    results = []
    for r, impl in schedule:
        res = run_single_trial(args.case_id, impl, r, r_script, run_dir, extra)
        results.append(res)
        if res["rc"] != 0:
            raise RuntimeError(f"Failed trial {res}")

    summary = {"case_id": args.case_id, "trials": results}
    for impl in ("native", "candidate"):
        sub = [x for x in results if x["impl"] == impl]
        summary[f"{impl}_wall_median"] = statistics.median(x["wall_seconds"] for x in sub)
        summary[f"{impl}_pss_median"] = statistics.median(x["peak_pss_gib"] for x in sub)
    summary["speedup_median_ratio"] = summary["native_wall_median"] / max(summary["candidate_wall_median"], 1e-9)
    summary["pss_change_pct"] = (summary["candidate_pss_median"] / max(summary["native_pss_median"], 1e-9) - 1.0) * 100.0
    (run_dir / f"{args.case_id}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
