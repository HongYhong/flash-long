"""Node 5 (T16-T17) tests: sparse Cython histograms, L1-cache float32 GBDT inference, CAGE bisect, PolyAFinder, TerminalCounters."""

from __future__ import annotations

import os
import sys
import types
from argparse import Namespace
from collections import Counter, namedtuple
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from isoquant_next.assignment.types import ReadAssignmentType
from isoquant_next.terminal import (
    EMPTY_COLUMNS,
    EMPTY_COLUMNS_GROUPED,
    HISTOGRAM_PAD,
    POLYA_MODEL_PATH,
    TSS_MODEL_PATH,
    CagePeakFinder,
    L1TreeEnsemble,
    PolyACounter,
    PolyAFinder,
    TSSCounter,
    TerminalHistogram,
    detect_peaks,
    detect_peaks_batch,
    get_polya_model,
    get_tss_model,
    rank_peaks_summary,
    reference_eval_json_tree_py,
    window_count,
)
from isoquant_next.terminal._cy_hist_stats import sparse_hist_moments_and_pad_cy

REPO_ROOT = Path(__file__).resolve().parents[2]
TESTS_DATA = REPO_ROOT / "isoquant_tests"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class _StubModel:
    def __init__(self, accept: bool = True) -> None:
        self._accept = accept

    def load_model(self, _path):
        return self

    def predict(self, X):
        value = 1 if self._accept else 0
        return np.full(len(X), value, dtype=np.int32)


@pytest.fixture
def stub_model(monkeypatch):
    import isoquant_next.terminal.counters as tc_mod

    m = _StubModel(accept=True)
    monkeypatch.setattr(tc_mod, "L1TreeEnsemble", lambda *a, **kw: m)
    return m


def _make_args(**extra):
    base = dict(genedb=True, fl_data=True, read_group=None)
    base.update(extra)
    return Namespace(**base)


def _make_read_assignment(
    *,
    strand: str = "+",
    chr_id: str = "chr1",
    polya_pos: int = 420,
    exons=None,
    assignment_type=ReadAssignmentType.unique,
    polyA_found: bool = True,
    transcript_id: str = "T1",
    gene_id: str = "G1",
    read_group_ids=None,
    gene_info=None,
):
    exons = exons or [(100, 200), (300, polya_pos)]
    polya_info = types.SimpleNamespace(
        external_polya_pos=polya_pos if strand == "+" else -1,
        external_polyt_pos=polya_pos if strand == "-" else -1,
        internal_polya_pos=-1,
        internal_polyt_pos=-1,
    )
    return types.SimpleNamespace(
        assignment_type=assignment_type,
        strand=strand,
        chr_id=chr_id,
        polyA_found=polyA_found,
        polya_info=polya_info,
        corrected_exons=exons,
        isoform_matches=[types.SimpleNamespace(assigned_transcript=transcript_id, assigned_gene=gene_id)],
        read_group_ids=list(read_group_ids) if read_group_ids else [],
        gene_info=gene_info or types.SimpleNamespace(all_isoforms_exons={transcript_id: list(exons)}),
    )


class TestCythonHistAndL1TreeEngine:
    def test_sparse_hist_moments_matches_scipy_and_numpy(self):
        rng = np.random.default_rng(42)
        test_vectors = [
            {100000100: 1},
            {100000100: 5, 100000105: 5},
            {420: 12, 421: 18, 422: 18, 500: 7, 503: 1},
        ]
        for _ in range(25):
            base = int(rng.integers(100, 200_000_000))
            offsets = rng.integers(0, 80, size=int(rng.integers(2, 20)))
            cnts = rng.integers(1, 50, size=len(offsets))
            d = {}
            for off, c in zip(offsets, cnts):
                d[base + int(off)] = d.get(base + int(off), 0) + int(c)
            test_vectors.append(d)

        for pos_counts in test_vectors:
            data = []
            for p, c in pos_counts.items():
                data.extend([p] * c)
            data_min = min(data)
            data_max = max(data)
            ref_hist = np.histogram(data, bins=(data_max - data_min) + 1, range=(data_min, data_max))[0]
            ref_padded = np.asarray([0] * HISTOGRAM_PAD + list(ref_hist) + [0] * HISTOGRAM_PAD, dtype=np.int32)
            ref_var = float(np.var(data))
            ref_mode_off = int(stats.mode(data, keepdims=False).mode) - data_min
            ref_skew = float(stats.skew(data)) if len(data) > 3 and ref_var >= 1e-12 else float("nan")
            ref_ent = float(stats.entropy(ref_padded))
            ref_mh = float(np.mean(ref_padded))

            items = sorted(pos_counts.items())
            u_pos = np.asarray([p for p, _ in items], dtype=np.int64)
            u_cnt = np.asarray([c for _, c in items], dtype=np.int64)
            cy_min, cy_n, cy_mode, cy_var, cy_skew, cy_ent, cy_mh, cy_pad = sparse_hist_moments_and_pad_cy(
                u_pos, u_cnt, HISTOGRAM_PAD
            )
            assert cy_min == data_min
            assert cy_n == len(data)
            assert cy_mode == ref_mode_off
            assert cy_var == pytest.approx(ref_var, rel=1e-9, abs=1e-9)
            if np.isnan(ref_skew):
                assert np.isnan(cy_skew)
            else:
                assert cy_skew == pytest.approx(ref_skew, rel=1e-7, abs=1e-7)
            assert cy_ent == pytest.approx(ref_ent, rel=1e-9, abs=1e-9)
            assert cy_mh == pytest.approx(ref_mh, rel=1e-12, abs=1e-12)
            assert np.array_equal(cy_pad, ref_padded)

    @pytest.mark.parametrize("model_path,get_fn", [(POLYA_MODEL_PATH, get_polya_model), (TSS_MODEL_PATH, get_tss_model)])
    def test_l1_tree_eval_matches_reference_json_under_1e_6(self, model_path, get_fn):
        ens = get_fn()
        assert isinstance(ens, L1TreeEnsemble)
        assert ens.num_raw_trees == 417
        assert 0 < ens.num_active_trees < 150  # Pruned single-leaf trees folded into const_margin

        rng = np.random.default_rng(123)
        X = np.column_stack(
            [
                rng.exponential(50.0, size=300),  # var
                rng.normal(0.0, 2.0, size=300),  # skew
                rng.integers(1, 15, size=300).astype(float),  # peak_count
                rng.uniform(0.5, 40.0, size=300),  # peak_width
                rng.uniform(0.0, 4.0, size=300),  # entropy
                rng.uniform(0.01, 20.0, size=300),  # mean_height
                rng.integers(1, 200, size=300).astype(float),  # peak_heights
                rng.uniform(0.05, 1.0, size=300),  # relative_height
            ]
        ).astype(np.float32)
        # Include NaN skew rows to verify default_left NaN routing
        X[:25, 1] = np.nan

        # Also construct real peak feature vectors (tight peak + noisy multi-peak)
        tight = TerminalHistogram({420: 15, 421: 30, 422: 12}).summarize()
        assert tight is not None and tight.peak_count >= 1
        rk = rank_peaks_summary(tight)[0]
        tight_row = np.asarray(
            [[tight.var, tight.skew, tight.peak_count, rk["peak_width"], tight.entropy, tight.mean_height, rk["peak_heights"], rk["relative_height"]]],
            dtype=np.float32,
        )
        X = np.vstack([X, tight_row])

        p_cy = ens.predict_proba_positive(X)
        p_ref = reference_eval_json_tree_py(model_path, X)
        max_diff = float(np.max(np.abs(p_cy - p_ref)))
        assert max_diff < 1e-6
        assert np.array_equal(ens.predict(X), (p_ref > 0.5).astype(np.int32))
        assert int(ens.predict(tight_row)[0]) == 1

        xgb_mod = pytest.importorskip("xgboost")
        clf = xgb_mod.XGBClassifier(n_jobs=1)
        clf.load_model(str(model_path))
        p_xgb = clf.predict_proba(X)[:, 1]
        assert float(np.max(np.abs(p_cy - p_xgb))) < 1e-6
        assert np.array_equal(ens.predict(X), clf.predict(X).astype(np.int32))


class TestPolyAAndCageOracle:
    def test_polya_finder_matches_upstream_oracle(self):
        from isoquant_lib.terminal_prediction.polya_finder import PolyAFinder as LibPF

        PolyAAlignment = namedtuple("PolyAAlignment", ("query_name", "cigartuples", "seq", "reference_start", "reference_end"))
        TUPLES1 = [(4, 2), (0, 15), (1, 1), (0, 5), (2, 1), (0, 16), (3, 4642), (0, 5), (1, 1), (0, 1), (1, 1), (0, 15), (2, 1), (0, 4), (1, 1), (0, 6), (2, 1), (0, 19), (1, 1), (0, 32), (1, 1), (0, 3), (3, 7649), (0, 25), (1, 1), (0, 13), (1, 1), (0, 29)]
        TUPLES2 = [(0, 167), (4, 30)]
        TUPLES3 = [(0, 167), (4, 30), (5, 5)]
        TUPLES4 = [(0, 150), (3, 100), (0, 10), (3, 100), (0, 11), (4, 27)]
        TUPLES5 = [(0, 147), (3, 100), (0, 5), (2, 2), (0, 4), (3, 100), (0, 5), (1, 2), (0, 5), (4, 30)]
        TUPLES6 = [(0, 147), (3, 100), (0, 5), (2, 2), (0, 4), (3, 100), (0, 5), (1, 2), (0, 5), (4, 30), (5, 5)]
        SEQ1 = "CTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGATGGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCGAAGAAGCAGGCATCGGAGACACTCCCAAACCAGGAGCGACCAAGCACTGGCGCATGTGACTCAAG"
        SEQ2 = "CTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGATGGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCGAAGAAGCAGGCATCGGACDCCTTTTTTTTTTAGGCGCCAAAAAAAAAAAACCAAAAAAAAAAAAA"
        SEQ3 = "TTTTTTTTTTTTTTTTTTTTCTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGATGGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCGAAGAAGCAGGCATCGGACGCCAAAAAAAAACCAAAAAAAAAAAAA"
        SEQ4 = "TTTTTTTATTTGTTTTTTTTTTTTTTTCTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGATGGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCGAAGAAGCAGGCATCGGACGCCAAAAAAAAACCAAAAAAAAAAAAA"

        lib_pf = LibPF()
        nxt_pf = PolyAFinder()
        alns = [
            PolyAAlignment("a1", TUPLES1, SEQ1, 51054, 63535),
            PolyAAlignment("a2", TUPLES2, SEQ2, 51054, 63535),
            PolyAAlignment("a3", TUPLES3, SEQ2, 51054, 63535),
            PolyAAlignment("a4", TUPLES4, SEQ2, 1000, 1300),
            PolyAAlignment("a5", TUPLES5, SEQ2, 1000, 1300),
            PolyAAlignment("a6", list(reversed(TUPLES2)), SEQ3, 51054, 63535),
            PolyAAlignment("a7", list(reversed(TUPLES3)), SEQ3, 51054, 63535),
            PolyAAlignment("a8", list(reversed(TUPLES5)), SEQ4, 1000, 1300),
            PolyAAlignment("a9", list(reversed(TUPLES6)), SEQ4, 1000, 1300),
        ]
        for aln in alns:
            p_lib = lib_pf.detect_polya(aln)
            p_nxt = nxt_pf.detect_polya(aln)
            assert (p_nxt.external_polya_pos, p_nxt.external_polyt_pos, p_nxt.internal_polya_pos, p_nxt.internal_polyt_pos) == (
                p_lib.external_polya_pos,
                p_lib.external_polyt_pos,
                p_lib.internal_polya_pos,
                p_lib.internal_polyt_pos,
            )

    def test_pure_python_cage_peak_finder_matches_pybedtools_oracle(self):
        from isoquant_lib.terminal_prediction.cage_finder import CagePeakFinder as LibCPF

        CageAlignment = namedtuple("CageAlignment", ("query_name", "reference_name", "reference_start", "reference_end", "is_reverse"))
        bed_path = TESTS_DATA / "toy_data" / "MAPT.Mouse.CAGE.bed"
        lib_cpf = LibCPF(str(bed_path), 50)
        nxt_cpf = CagePeakFinder(bed_path, 50)

        queries = [
            CageAlignment("q1", "chr11", 24, 2165, False),
            CageAlignment("q2", "chr11", 0, 1296, False),
            CageAlignment("q3", "chr11", 80, 372, True),
            CageAlignment("q4", "chr11", 57, 522, False),
            CageAlignment("q5", "chr11", 69, 1499, False),
            CageAlignment("q6", "chr11", 134, 1322, False),
            CageAlignment("q7", "chr11", 102, 2191, False),
            CageAlignment("q8", "chr11", 11, 126, True),
        ]
        for q in queries:
            lib_hits = [(h.chrom, h.start, h.end, h.name) for h in lib_cpf.find_cage_peak(q)]
            nxt_hits = [(h.chrom, h.start, h.end, h.name) for h in nxt_cpf.find_cage_peak(q)]
            assert sorted(nxt_hits) == sorted(lib_hits)


class TestTerminalPeaksAndCounters:
    @pytest.mark.parametrize(
        "positions",
        [
            [420],
            [420 + i % 3 for i in range(30)],
            [420 + i % 3 for i in range(30)] + [500 + i % 3 for i in range(20)],
            [400, 401],
        ],
    )
    def test_detect_peaks_parity_with_counter(self, stub_model, tmp_path, positions):
        out = tmp_path / "parity.tsv"
        counter = PolyACounter(_make_args(), str(out))
        exons = [(100, 200), (300, 400)]
        for pos in positions:
            counter.add_read_info(_make_read_assignment(polya_pos=pos, exons=exons))
        counter.dump()
        df = pd.read_csv(out, sep="\t")
        counter_peaks = set() if df.empty else set(zip(df["prediction"].astype(int), df["counts"].astype(int)))
        det = {(p.position, p.count) for p in detect_peaks(dict(Counter(positions)), _StubModel(accept=True))}
        assert det == counter_peaks

    def test_detect_peaks_empty_batch_and_rejected(self):
        assert detect_peaks({}, _StubModel(accept=True)) == []
        hists = [{420: 1}, {420: 10, 421: 10, 422: 10}, {}]
        m = _StubModel(accept=True)
        assert detect_peaks_batch(hists, m) == [detect_peaks(h, m) for h in hists]
        assert {(p.position, p.count) for p in detect_peaks(Counter([420 + i % 3 for i in range(30)]), _StubModel(accept=False))} == set()

    def test_grouped_searchsorted_and_flush_invariants(self, stub_model, tmp_path):
        pool = types.SimpleNamespace(get_str=lambda i: f"g{i}")
        string_pools = types.SimpleNamespace(get_read_group_pool=lambda _idx: pool)
        out = tmp_path / "grouped.tsv"
        counter = PolyACounter(_make_args(), str(out), string_pools=string_pools, group_index=0)
        exons = [(100, 200), (300, 400)]
        for offset in range(30):
            counter.add_read_info(_make_read_assignment(polya_pos=420 + (offset % 3), exons=exons, read_group_ids=[offset % 2]))
        counter.flush()
        assert counter.transcripts
        assert counter.last_gene_predictions == []
        counter.dump()
        df = pd.read_csv(out, sep="\t")
        assert list(df.columns) == EMPTY_COLUMNS_GROUPED
        assert set(df["group_id"]) == {"g0", "g1"}
        assert int(df["counts_byGroup"].sum()) == 30

    def test_per_gene_flush_and_known_novel_and_rejection(self, monkeypatch, tmp_path):
        import isoquant_next.terminal.counters as tc_mod

        monkeypatch.setattr(tc_mod, "L1TreeEnsemble", lambda *a, **kw: _StubModel(accept=True))
        out = tmp_path / "split.tsv"
        counter = PolyACounter(_make_args(), str(out))
        exons = [(100, 200), (300, 400)]
        for offset in range(10):
            counter.add_read_info(_make_read_assignment(polya_pos=420 + offset % 3, exons=exons))
        counter.flush()
        assert counter.last_gene_predictions
        assert counter.transcripts == {}
        assert counter._all_transcripts
        for offset in range(20):
            counter.add_read_info(_make_read_assignment(polya_pos=420 + offset % 3, exons=exons))
        counter.flush()
        counter.dump()
        df = pd.read_csv(out, sep="\t")
        assert list(df.columns) == EMPTY_COLUMNS
        assert int(df["counts"].sum()) == 30
        assert (df["flag"] == "Novel").all()

        # Known within tolerance
        known_out = tmp_path / "known.tsv"
        c_known = PolyACounter(_make_args(), str(known_out))
        for _ in range(30):
            c_known.add_read_info(_make_read_assignment(polya_pos=402, exons=exons))
        c_known.dump()
        assert (pd.read_csv(known_out, sep="\t")["flag"] == "Known").all()

        # TSS counter strand coordinates
        t_out = tmp_path / "tss.tsv"
        c_tss = TSSCounter(_make_args(), str(t_out))
        c_tss.add_read_info(_make_read_assignment(strand="+", exons=[(105, 200), (300, 400)]))
        c_tss.add_read_info(_make_read_assignment(strand="-", exons=[(100, 200), (300, 405)]))
        assert sorted(c_tss.transcripts["T1"].histogram.counts.keys()) == [105, 405]
        c_tss.add_read_info(_make_read_assignment(strand="."))
        assert sorted(c_tss.transcripts["T1"].histogram.counts.keys()) == [105, 405]

        # Model rejection drops peak
        monkeypatch.setattr(tc_mod, "L1TreeEnsemble", lambda *a, **kw: _StubModel(accept=False))
        rej_out = tmp_path / "rej.tsv"
        c_rej = PolyACounter(_make_args(), str(rej_out))
        for offset in range(30):
            c_rej.add_read_info(_make_read_assignment(polya_pos=420 + offset % 3, exons=exons))
        c_rej.dump()
        assert pd.read_csv(rej_out, sep="\t").empty

    def test_polya_counter_tsv_exact_parity_with_oracle(self, monkeypatch, tmp_path):
        import isoquant_next.terminal.counters as nxt_tc

        monkeypatch.setitem(sys.modules, "xgboost", types.SimpleNamespace(XGBClassifier=_StubModel))
        from isoquant_lib.assignment.isoform_assignment import ReadAssignmentType as LibRAT
        from isoquant_lib.terminal_prediction import terminal_counter as lib_tc

        monkeypatch.setattr(nxt_tc, "L1TreeEnsemble", lambda *a, **kw: _StubModel(accept=True))
        monkeypatch.setattr(lib_tc, "XGBClassifier", lambda *a, **kw: _StubModel(accept=True))

        exons = [(100, 200), (300, 400)]
        lib_out = tmp_path / "lib.tsv"
        nxt_out = tmp_path / "nxt.tsv"
        pool = types.SimpleNamespace(get_str=lambda i: f"g{i}")
        sp = types.SimpleNamespace(get_read_group_pool=lambda _idx: pool)

        c_lib = lib_tc.PolyACounter(_make_args(), str(lib_out), string_pools=sp, group_index=0)
        c_nxt = PolyACounter(_make_args(), str(nxt_out), string_pools=sp, group_index=0)
        for offset in range(30):
            ra_nxt = _make_read_assignment(polya_pos=420 + (offset % 3), exons=exons, read_group_ids=[offset % 2])
            ra_lib = _make_read_assignment(
                polya_pos=420 + (offset % 3),
                exons=exons,
                read_group_ids=[offset % 2],
                assignment_type=LibRAT.unique,
            )
            c_lib.add_read_info(ra_lib)
            c_nxt.add_read_info(ra_nxt)
        c_lib.dump()
        c_nxt.dump()
        assert nxt_out.read_text() == lib_out.read_text()
