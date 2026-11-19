from __future__ import annotations

import os
import subprocess
import sys
from collections import namedtuple
from enum import Enum
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from isoquant_next.config import (
    AnalysisType,
    IsoQuantMode,
    PipelineConfig,
)
from isoquant_next.core.concurrency import WorkerPoolExecutor
from isoquant_next.core.error_codes import IsoQuantError, IsoQuantExitCode
from isoquant_next.core.file_naming import (
    convert_chr_id_to_file_name_str,
    counts_file_name,
    filtered_reads_file_name,
    multimappers_file_name,
    saves_file_name,
)
from isoquant_next.core.file_utils import merge_file_list, merge_files, normalize_path
from isoquant_next.core.id_policy import (
    AtomicIDDistributor,
    ExcludingIdDistributor,
    SimpleIDDistributor,
)
from isoquant_next.core.logging_setup import pin_single_thread_math_env
from isoquant_next.core.stats import EnumStats
from isoquant_next.core.utils import (
    get_first_best_from_sorted,
    get_path_to_program,
    interval_bin_search,
    interval_bin_search_rev,
    interval_len,
    intervals_total_length,
    jaccard_similarity,
    proper_plural_form,
    reverse_complement,
    rindex,
    rreplace,
    sum_intervals_from_point,
    sum_intervals_to_point,
)


def _worker_square(x: int) -> tuple[int, str | None]:
    return x * x, os.environ.get("OMP_NUM_THREADS")


def test_pipeline_config_resolve_and_manifest_roundtrip(tmp_path: Path) -> None:
    cfg = PipelineConfig.resolve(
        data_type="ont",
        mode=IsoQuantMode.tenX_v3,
        genedb="/tmp/fake.gtf",
        threads=4,
        output_dir=str(tmp_path / "out"),
    )
    assert cfg.data_type == "nanopore"
    assert cfg.matching.strategy_name == "default"
    assert cfg.matching.delta == 6
    assert cfg.barcodes.read_group == ("barcode",)
    assert cfg.analyses == (AnalysisType.quantification,)

    manifest = tmp_path / "run_manifest.json"
    cfg.save_manifest(manifest)
    loaded = PipelineConfig.load_manifest(manifest)
    assert loaded == cfg
    assert loaded.command_fingerprint() == cfg.command_fingerprint()


def test_pipeline_config_rejects_negative_delta() -> None:
    with pytest.raises(IsoQuantError) as exc_info:
        PipelineConfig.resolve(data_type="pacbio", delta=-1)
    assert exc_info.value.code == IsoQuantExitCode.INVALID_PARAMETER


def test_reverse_complement_supports_iupac_and_case() -> None:
    assert reverse_complement("AcGtUuNnRYSWKM") == "KMWSRYnNaAaCgT"


def test_get_path_to_program_handles_major_version_upgrade(monkeypatch) -> None:
    def _mock_run(*args, **kwargs):
        cp = namedtuple("CompletedProcess", ("stdout",))
        return cp(b"tool 2.1.0\n")

    monkeypatch.setattr(subprocess, "run", _mock_run)
    monkeypatch.setattr(os.path, "isfile", lambda p: True)
    monkeypatch.setattr(os, "access", lambda p, m: True)

    assert get_path_to_program("my_prog", "/folder", "1.38") == "/folder/my_prog"
    assert get_path_to_program("my_prog", "/folder", "2.2") is None


def test_interval_binary_search_exact_and_gaps() -> None:
    intervals = [(10, 20), (30, 40), (50, 60)]
    assert interval_bin_search(intervals, 35) == 1
    assert interval_bin_search(intervals, 25) == 0
    assert interval_bin_search_rev(intervals, 25) == 1
    assert interval_bin_search_rev(intervals, 15) == 0
    assert interval_bin_search(intervals, 5) == -1
    assert interval_bin_search(intervals, 100) == -1
    assert interval_bin_search([], 100) == -1
    assert interval_bin_search_rev([], 100) == -1


def test_range_math_and_utilities() -> None:
    assert proper_plural_form("read", 1) == "1 read"
    assert proper_plural_form("read", 2) == "2 reads"
    assert rreplace("a_b_b", "b", "c") == "a_b_c"
    assert rindex([1, 1, 0], 1) == 1
    assert get_first_best_from_sorted([(1, 1), (2, 1), (3, 3)]) == [1, 2]
    assert interval_len((1, 20)) == 20
    assert intervals_total_length([(1, 20), (40, 50)]) == 31
    assert sum_intervals_to_point([(1, 20), (40, 50)], 41) == 21
    assert sum_intervals_from_point([(1, 20), (40, 50)], 39) == 11
    assert pytest.approx(jaccard_similarity([(1, 70)], [(31, 100)])) == 0.4


def test_file_naming_and_merge_collision(tmp_path: Path) -> None:
    assert "chr1_2" in convert_chr_id_to_file_name_str("chr1/2")
    assert "multimappers" in multimappers_file_name("out", "chr1")
    assert "filtered" in filtered_reads_file_name("out", "chr1:10-20")
    assert saves_file_name("out", "chr1") == "out_chr1"
    assert counts_file_name("p", linear=True) == "p.linear.tsv"

    assert merge_file_list("/out/p/p.gene_counts.tsv", "p", ["chr9"]) == [
        "/out/p/p_chr9.gene_counts.tsv"
    ]
    assert normalize_path("/a/b/cfg.yaml", "rel/x.bam") == "/a/b/rel/x.bam"

    for ch, val in [("chr1", "a\n"), ("chr2", "b\n")]:
        (tmp_path / f"f_{ch}.tsv").write_text(f"#hdr\n{val}", encoding="utf-8")
    merged = tmp_path / "f.tsv"
    with open(merged, "w", encoding="utf-8") as fh:
        merge_files(str(merged), "f", ["chr1", "chr2"], fh, copy_header=True)
    assert merged.read_text(encoding="utf-8") == "#hdr\na\nb\n"


def test_id_distributors() -> None:
    s = SimpleIDDistributor()
    assert [s.increment() for _ in range(3)] == [1, 2, 3]

    mock_genedb = MagicMock()
    g = MagicMock()
    g.id = "novel_gene_5"
    t = MagicMock()
    t.id = "transcript10.1"
    mock_genedb.region = lambda seqid, start, featuretype: [g] if featuretype == "gene" else [t]
    excl = ExcludingIdDistributor(mock_genedb, "chr1")
    ids = [excl.increment() for _ in range(10)]
    assert 5 not in ids and 10 not in ids
    assert ids == [1, 2, 3, 4, 6, 7, 8, 9, 11, 12]

    atom = AtomicIDDistributor()
    assert atom.increment() == 1


def test_enum_stats_roundtrip(tmp_path: Path) -> None:
    class Demo(Enum):
        A = 1
        B = 2

    st = EnumStats()
    st.add(Demo.A, 3)
    st.add(Demo.B, 5)
    out = tmp_path / "stats.tsv"
    st.dump(str(out))
    reloaded = EnumStats()
    reloaded.load(str(out), enum_cls=Demo)
    assert reloaded.stats_dict[Demo.A] == 3
    assert reloaded.stats_dict[Demo.B] == 5


def test_worker_pool_executor_pins_omp_and_runs_sequential_and_parallel() -> None:
    pin_single_thread_math_env(force=True)
    seq_res = WorkerPoolExecutor(max_workers=1).map_ordered(_worker_square, [2, 3])
    assert seq_res == [(4, "1"), (9, "1")]

    par_res = WorkerPoolExecutor(max_workers=2, mp_method="spawn").map_ordered(_worker_square, [2, 3, 4])
    assert par_res == [(4, "1"), (9, "1"), (16, "1")]
