"""Node 10 (T23-T24) tests: Visualization (`OutputConfig`, `.iqidx`, Polars, Agg) & NumPy 2/Pandas 3 compat."""

from __future__ import annotations

import ast
import json
import os
from pathlib import Path
import pickle
import sys
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

from isoquant_next.annotation.iqidx import compile_gtf_to_iqidx
from isoquant_next.barcodes.common import batch_str_to_2bit, bit_to_str, str_to_2bit
from isoquant_next.barcodes.indexers import Array2BitKmerIndexer
from isoquant_next.config.settings import PipelineConfig
from isoquant_next.visualization.cli import run_visualize
from isoquant_next.visualization.dict_builder import DictionaryBuilder
from isoquant_next.visualization.output_config import OutputConfig
from isoquant_next.visualization.plotter import PlotOutput

REFACTOR_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = REFACTOR_ROOT.parent
MAPT_GTF = REPO_ROOT / "isoquant_tests" / "toy_data" / "MAPT.Mouse.genedb.gtf"
SRC_PKG_DIR = REFACTOR_ROOT / "src" / "isoquant_next"


class _LegacyParamsStub:
    def __init__(self, genedb: str, prefix: str = "LEG", data_type: str = "nanopore") -> None:
        self.genedb = genedb
        self.original_annotation = genedb
        self.prefix = prefix
        self.data_type = data_type


def _current_rss_kib() -> int:
    with open("/proc/self/status", "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("VmRSS:"):
                return int(line.split()[1])
    return 0


class TestVisualizationArchitecture:
    def test_output_config_manifest_priority_and_legacy_params_fallback(self, tmp_path: Path):
        # 1. Both run_manifest.json and .params exist -> run_manifest.json takes precedence
        out_dir = tmp_path / "manifest_run"
        sample_dir = out_dir / "S_MAN"
        sample_dir.mkdir(parents=True)
        iqidx_path = out_dir / "MAPT.Mouse.genedb.iqidx"
        compile_gtf_to_iqidx(MAPT_GTF, iqidx_path)

        cfg = PipelineConfig.resolve(
            data_type="pacbio_ccs",
            genedb=str(iqidx_path),
            output_dir=str(out_dir),
            prefix="S_MAN",
        )
        cfg.save_manifest(out_dir / "run_manifest.json")
        (out_dir / ".params").write_bytes(pickle.dumps(_LegacyParamsStub(str(MAPT_GTF), prefix="WRONG", data_type="nanopore")))
        (sample_dir / "S_MAN.gene_counts.tsv").write_text("#feature_id\tcount\nG1\t10\n", encoding="utf-8")
        (sample_dir / "S_MAN.transcript_models.gtf").write_text("", encoding="utf-8")

        loaded = OutputConfig(str(out_dir))
        assert loaded.source_kind == "run_manifest"
        assert loaded.data_type == "pacbio_ccs"
        assert loaded.iqidx_path == iqidx_path
        assert len(loaded.samples) == 1
        assert loaded.samples[0].prefix == "S_MAN"
        assert loaded.ref_only is False

        # 2. Only legacy .params exists -> read-only fallback succeeds without crashing
        leg_dir = tmp_path / "legacy_run"
        leg_sample = leg_dir / "LEG"
        leg_sample.mkdir(parents=True)
        (leg_dir / ".params").write_bytes(pickle.dumps(_LegacyParamsStub(str(MAPT_GTF), prefix="LEG", data_type="nanopore")))
        (leg_sample / "LEG.gene_counts.tsv").write_text("#feature_id\tcount\nG1\t5\n", encoding="utf-8")

        legacy_cfg = OutputConfig(str(leg_dir))
        assert legacy_cfg.source_kind == "legacy_params"
        assert legacy_cfg.reference_gtf == MAPT_GTF
        assert legacy_cfg.samples[0].prefix == "LEG"
        assert legacy_cfg.ref_only is True

    def test_dictionary_builder_iqidx_and_streaming_gtf_without_gffutils(self, tmp_path: Path):
        out_dir = tmp_path / "viz_dict"
        sdir = out_dir / "MAPT"
        sdir.mkdir(parents=True)
        iqidx_path = out_dir / "MAPT.iqidx"
        compile_gtf_to_iqidx(MAPT_GTF, iqidx_path)

        (sdir / "MAPT.transcript_models.gtf").write_text(
            'chr11\tIsoQuant\texon\t50979\t51091\t.\t+\t.\tgene_id "ENSMUSG00000018411.17"; transcript_id "novel_tx_1";\n'
            'chr11\tIsoQuant\texon\t63470\t63535\t.\t+\t.\tgene_id "ENSMUSG00000018411.17"; transcript_id "novel_tx_1";\n',
            encoding="utf-8",
        )
        cfg = OutputConfig(str(out_dir), gtf=str(iqidx_path))
        builder = DictionaryBuilder(cfg)
        genes_iqidx = builder.build()
        assert "ENSMUSG00000018411.17" in genes_iqidx
        assert "novel_tx_1" in genes_iqidx["ENSMUSG00000018411.17"].transcripts
        assert genes_iqidx["ENSMUSG00000018411.17"].transcripts["novel_tx_1"].is_discovered is True

        # Streaming GTF path without compiling .db or importing gffutils
        cfg_gtf = OutputConfig(str(out_dir), gtf=str(MAPT_GTF), ref_only=True)
        cfg_gtf.iqidx_path = None
        genes_gtf = DictionaryBuilder(cfg_gtf).build()
        assert set(genes_gtf.keys()) == set(genes_iqidx.keys())
        assert "novel_tx_1" not in genes_gtf["ENSMUSG00000018411.17"].transcripts

    def test_plotter_agg_backend_close_and_memory_stability(self, tmp_path: Path):
        assert matplotlib.get_backend().lower() == "agg"
        out_dir = tmp_path / "plots"
        iqidx_path = tmp_path / "MAPT.iqidx"
        compile_gtf_to_iqidx(MAPT_GTF, iqidx_path)
        (tmp_path / "S").mkdir()
        cfg = OutputConfig(str(tmp_path), gtf=str(iqidx_path), ref_only=True)
        gmap = DictionaryBuilder(cfg).build()
        rec = gmap["ENSMUSG00000018411.17"]
        counts = {tid: {"G1": float((i + 1) * 3)} for i, tid in enumerate(rec.transcripts)}

        plotter = PlotOutput(out_dir, image_format="png", dpi=72)
        p0 = plotter.plot_gene_model(rec, counts_by_tx=counts)
        assert p0 is not None and p0.exists() and p0.stat().st_size > 500
        p_pie = plotter.plot_count_summary_pie({"unique": 40, "ambiguous": 10}, filename="stats")
        assert p_pie is not None and p_pie.exists()
        assert plt.get_fignums() == []

        # Warm up FreeType/matplotlib glyph caches first, then assert steady-state RSS delta across 50 renders
        for warmup_idx in range(5):
            rec.gene_id = f"ENSMUSG_WARM_{warmup_idx}"
            plotter.plot_gene_model(rec, counts_by_tx=counts)
            assert plt.get_fignums() == []

        rss_before = _current_rss_kib()
        for idx in range(50):
            rec.gene_id = f"ENSMUSG_LOOP_{idx % 5}"
            plotter.plot_gene_model(rec, counts_by_tx=counts)
            assert plt.get_fignums() == []
        rss_after = _current_rss_kib()
        assert (rss_after - rss_before) < 10 * 1024

    def test_gene_ranking_polars_and_visualize_cli_without_pandas(self, tmp_path: Path):
        pytest.importorskip("polars")
        from isoquant_next.visualization.gene_ranking import rank_and_visualize_genes

        out_dir = tmp_path / "viz_cli"
        sdir = out_dir / "S1"
        sdir.mkdir(parents=True)
        iqidx_path = out_dir / "MAPT.iqidx"
        compile_gtf_to_iqidx(MAPT_GTF, iqidx_path)

        cfg = OutputConfig(str(out_dir), gtf=str(iqidx_path), ref_only=True)
        gmap = DictionaryBuilder(cfg).build()
        tx_ids = list(gmap["ENSMUSG00000018411.17"].transcripts.keys())
        assert len(tx_ids) >= 2
        tx_to_gene = {tid: gid for gid, g in gmap.items() for tid in g.transcripts}

        counts_tsv = sdir / "S1.transcript_counts.tsv"
        lines = ["#feature_id\tcond_A\tcond_B"]
        lines.append(f"{tx_ids[0]}\t100.0\t5.0")
        lines.append(f"{tx_ids[1]}\t5.0\t100.0")
        for extra in tx_ids[2:]:
            lines.append(f"{extra}\t1.0\t1.0")
        counts_tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")

        ranked = rank_and_visualize_genes(counts_tsv, tx_to_gene, num_genes=5, min_total_count=5.0)
        assert ranked and ranked[0].gene_id == "ENSMUSG00000018411.17"
        assert ranked[0].divergence_score > 0.5

        plot_out = tmp_path / "cli_plots"
        rc = run_visualize(
            [
                "-o",
                str(out_dir),
                "--gtf",
                str(iqidx_path),
                "--counts",
                "--find_genes",
                "--num_genes",
                "3",
                "--format",
                "png",
                "--visual_output",
                str(plot_out),
            ]
        )
        assert rc == 0
        assert any(plot_out.glob("*.png"))


class TestForwardCompatNumpy2AndPandas3:
    def test_zero_forbidden_runtime_imports_in_src(self):
        forbidden = {"xgboost", "ssw", "editdistance", "mappy", "pybedtools", "Bio", "pandas"}
        violations = []
        for py_file in sorted(SRC_PKG_DIR.rglob("*.py")):
            tree = ast.parse(py_file.read_text(encoding="utf-8"), filename=str(py_file))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        top = alias.name.split(".", 1)[0]
                        if top in forbidden:
                            violations.append(f"{py_file}: import {alias.name}")
                elif isinstance(node, ast.ImportFrom) and node.module:
                    top = node.module.split(".", 1)[0]
                    if top in forbidden:
                        violations.append(f"{py_file}: from {node.module}")
        for cy_file in sorted(list(SRC_PKG_DIR.rglob("*.pyx")) + list(SRC_PKG_DIR.rglob("*.pxd"))):
            for lineno, line in enumerate(cy_file.read_text(encoding="utf-8").splitlines(), start=1):
                stripped = line.strip()
                if stripped.startswith("#"):
                    continue
                for pkg in forbidden:
                    if f"import {pkg}" in stripped or f"from {pkg}" in stripped:
                        violations.append(f"{cy_file}:{lineno}: {stripped}")
        assert violations == []

    def test_numpy2_uint64_bitwise_and_kmer_compat(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            warnings.simplefilter("error", FutureWarning)

            k = 14
            mask = np.uint64((1 << (2 * k)) - 1)
            encoded = np.uint64(str_to_2bit("ACGTACGTACGTAC"))
            shifted = (encoded << np.uint64(2)) & mask
            assert shifted.dtype == np.uint64
            assert bit_to_str(int(encoded), k) == "ACGTACGTACGTAC"

            batch = batch_str_to_2bit(["AAACCCGGGTTTAA", "TTTAAACCCGGGTT"], seq_len=k)
            assert batch.dtype == np.uint64 and batch.shape == (2,)
            idx = Array2BitKmerIndexer(batch, kmer_size=6, seq_len=k)
            hits = idx.get_occurrences("AAACCCGGGTTTAA")
            assert len(hits) > 0 and hits[0][0] == "AAACCCGGGTTTAA"

    def test_pandas3_cow_and_scipy_mtx_forward_compat(self, tmp_path: Path):
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            warnings.simplefilter("error", FutureWarning)

            raw_df = pd.DataFrame(
                {
                    "feature_id": ["G1", "G2", "G3"],
                    "cell_A": [3.0, 0.0, 2.0],
                    "cell_B": [0.0, 4.0, 1.0],
                }
            )
            # Copy-on-Write safe transformation (explicit copy / new Series, no chained mutation)
            norm_df = raw_df.copy()
            mat_vals = norm_df[["cell_A", "cell_B"]].to_numpy(dtype=np.float64)
            col_totals = np.maximum(mat_vals.sum(axis=0, keepdims=True), 1.0)
            norm_vals = mat_vals / col_totals
            norm_df = norm_df.assign(cell_A=norm_vals[:, 0], cell_B=norm_vals[:, 1])
            assert abs(float(norm_df["cell_A"].sum()) - 1.0) < 1e-12
            assert abs(float(raw_df["cell_A"].sum()) - 5.0) < 1e-12

            rows, cols = np.nonzero(mat_vals > 0)
            coo = coo_matrix((mat_vals[rows, cols], (rows, cols)), shape=mat_vals.shape, dtype=np.float64)
            assert coo.nnz == 4
            assert coo.shape == (3, 2)
            mtx_path = tmp_path / "matrix.mtx"
            mmwrite(mtx_path, coo)
            assert mtx_path.exists() and mtx_path.stat().st_size > 0
