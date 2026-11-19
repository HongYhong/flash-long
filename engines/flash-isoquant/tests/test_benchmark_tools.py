from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts" / "benchmark"


def _load_script_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    import sys

    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_run_monitored_command_records_pss_and_rusage(tmp_path: Path) -> None:
    bench = _load_script_module("run_tiered_benchmark", SCRIPTS_DIR / "run_tiered_benchmark.py")
    metrics = bench.run_monitored_command(
        ["python", "-c", "import time, numpy as np; a = np.ones(2000000, dtype=np.float64); time.sleep(0.8); print(a.sum())"],
        label="smoke",
        rep=1,
        log_dir=tmp_path,
        interval=0.2,
    )
    assert metrics.exit_code == 0
    assert metrics.wall_seconds >= 0.7
    assert metrics.peak_tree_pss_kib > 0
    assert metrics.peak_tree_rss_kib >= metrics.peak_tree_pss_kib
    assert Path(metrics.pss_trace_path).exists()


def test_compare_outputs_counts_and_gtf(tmp_path: Path) -> None:
    cmp_mod = _load_script_module("compare_outputs", SCRIPTS_DIR / "compare_outputs.py")
    ref_tsv = tmp_path / "ref.tsv"
    cand_tsv = tmp_path / "cand.tsv"
    ref_tsv.write_text("#feature_id\ts1\nG1\t10\nG2\t2.500000\n", encoding="utf-8")
    cand_tsv.write_text("#feature_id\ts1\nG1\t10\nG2\t2.5000005\n", encoding="utf-8")
    res_tsv = cmp_mod.compare_count_tables(ref_tsv, cand_tsv, float_atol=1e-6)
    assert res_tsv.passed

    ref_gtf = tmp_path / "ref.gtf"
    cand_gtf = tmp_path / "cand.gtf"
    gtf_txt = (
        'chr1\tIsoQuant\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tIsoQuant\texon\t400\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    )
    ref_gtf.write_text(gtf_txt, encoding="utf-8")
    cand_gtf.write_text(gtf_txt, encoding="utf-8")
    res_gtf = cmp_mod.compare_gtf_models(ref_gtf, cand_gtf)
    assert res_gtf.passed
