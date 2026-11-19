"""Unified sparse TerminalHistogram & peak detector backed by Cython sparse stats & L1 GBDT."""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, NamedTuple, Optional, Tuple

import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths

from isoquant_next.terminal._cy_hist_stats import sparse_hist_moments_and_pad_cy
from isoquant_next.terminal.tree_model import FEATURE_COLUMNS, L1TreeEnsemble

logger = logging.getLogger("IsoQuant")

_DATA_DIR = Path(__file__).resolve().parent.parent / "data"
POLYA_MODEL_PATH = _DATA_DIR / "model_polya.json"
TSS_MODEL_PATH = _DATA_DIR / "model_tss.json"

PEAK_DISTANCE = 10
PEAK_REL_HEIGHT = 0.98
HISTOGRAM_PAD = 10
ANNOTATION_TOLERANCE = 10

_MODEL_CACHE: Dict[str, L1TreeEnsemble] = {}


def _get_model(model_path: Path) -> L1TreeEnsemble:
    key = str(model_path)
    model = _MODEL_CACHE.get(key)
    if model is None:
        model = L1TreeEnsemble(model_path)
        _MODEL_CACHE[key] = model
    return model


def get_polya_model() -> L1TreeEnsemble:
    return _get_model(POLYA_MODEL_PATH)


def get_tss_model() -> L1TreeEnsemble:
    return _get_model(TSS_MODEL_PATH)


class Peak(NamedTuple):
    """Accepted terminal-position peak in genomic coordinates."""

    position: int
    count: int
    left: int
    right: int


class RawHistogramSummary(NamedTuple):
    start: int
    mode: int
    var: float
    skew: float
    entropy: float
    mean_height: float
    padded: np.ndarray
    peak_count: int
    peak_location: List[int]
    peak_width: List[float]
    peak_prominence: List[float]
    peak_left: List[int]
    peak_right: List[int]
    peak_heights: List[int]


class TerminalHistogram:
    """Sparse `{genomic_pos: count}` accumulator eliminating `[pos] * count` list expansion."""

    __slots__ = ("counts", "_total")

    def __init__(self, initial: Optional[Mapping[int, int]] = None) -> None:
        self.counts: Dict[int, int] = defaultdict(int)
        self._total = 0
        if initial:
            for pos, cnt in initial.items():
                c = int(cnt)
                if c > 0:
                    self.counts[int(pos)] += c
                    self._total += c

    def add(self, pos: int, count: int = 1) -> None:
        if count > 0:
            p = int(pos)
            self.counts[p] += count
            self._total += count

    def extend(self, positions: Iterable[int]) -> None:
        for p in positions:
            self.counts[int(p)] += 1
            self._total += 1

    def merge(self, other: "TerminalHistogram") -> None:
        for pos, cnt in other.counts.items():
            if cnt > 0:
                self.counts[pos] += cnt
                self._total += cnt

    def __bool__(self) -> bool:
        return self._total > 0

    def __len__(self) -> int:
        return self._total

    def to_sorted_arrays(self) -> Tuple[np.ndarray, np.ndarray]:
        items = sorted((p, c) for p, c in self.counts.items() if c > 0)
        if not items:
            return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
        u_pos = np.fromiter((p for p, _ in items), dtype=np.int64, count=len(items))
        u_cnt = np.fromiter((c for _, c in items), dtype=np.int64, count=len(items))
        return u_pos, u_cnt

    def summarize(self, compute_prominence: bool = False) -> Optional[RawHistogramSummary]:
        u_pos, u_cnt = self.to_sorted_arrays()
        if u_pos.size == 0:
            return None
        data_min, _, mode_off, var, skew_val, entropy_val, mean_height, padded = sparse_hist_moments_and_pad_cy(
            u_pos, u_cnt, HISTOGRAM_PAD
        )
        peaks_idx = find_peaks(padded, distance=PEAK_DISTANCE)[0]
        peak_count = int(peaks_idx.shape[0])
        if peak_count == 0:
            return RawHistogramSummary(
                start=int(data_min),
                mode=int(mode_off),
                var=float(var),
                skew=float(skew_val),
                entropy=float(entropy_val),
                mean_height=float(mean_height),
                padded=padded,
                peak_count=0,
                peak_location=[],
                peak_width=[],
                peak_prominence=[],
                peak_left=[],
                peak_right=[],
                peak_heights=[],
            )

        prom = peak_prominences(padded, peaks_idx)[0].tolist() if compute_prominence else [0.0] * peak_count
        widths = peak_widths(padded, peaks_idx, rel_height=PEAK_REL_HEIGHT)
        peak_location = [int(j) - HISTOGRAM_PAD for j in peaks_idx]
        peak_left = [int(j) - HISTOGRAM_PAD for j in widths[2]]
        peak_right = [int(j) - HISTOGRAM_PAD for j in widths[3]]
        peak_heights = [int(padded[p + HISTOGRAM_PAD]) for p in peak_location]
        return RawHistogramSummary(
            start=int(data_min),
            mode=int(mode_off),
            var=float(var),
            skew=float(skew_val),
            entropy=float(entropy_val),
            mean_height=float(mean_height),
            padded=padded,
            peak_count=peak_count,
            peak_location=peak_location,
            peak_width=[float(w) for w in widths[0]],
            peak_prominence=[float(v) for v in prom],
            peak_left=peak_left,
            peak_right=peak_right,
            peak_heights=peak_heights,
        )


def rank_peaks_summary(summary: RawHistogramSummary) -> List[dict]:
    heights = np.asarray(summary.peak_heights, dtype=float)
    if heights.size > 1:
        order = np.argsort(-heights)
        top = float(heights[order[0]])
    else:
        order = [0]
        top = None

    ranked: List[dict] = []
    for r_i, idx in enumerate(order):
        h = summary.peak_heights[idx]
        ranked.append(
            {
                "peak_location": summary.peak_location[idx],
                "peak_width": summary.peak_width[idx],
                "peak_prominence": summary.peak_prominence[idx],
                "peak_left": summary.peak_left[idx],
                "peak_right": summary.peak_right[idx],
                "peak_heights": h,
                "rank": (r_i + 1) if heights.size > 1 else 0,
                "relative_height": (1.0 if top is None else (float(h) / top if top else 0.0)),
            }
        )
    return ranked


def window_count(padded: np.ndarray, left: int, right: int) -> int:
    low = max(0, int(left) + HISTOGRAM_PAD)
    high = min(int(padded.shape[0]), int(right) + HISTOGRAM_PAD + 1)
    if high <= low:
        return 0
    return int(padded[low:high].sum())


def detect_peaks(position_counts: Mapping[int, int], model=None) -> List[Peak]:
    """Detect accepted peaks in a single `{genomic_position: count}` histogram."""
    return detect_peaks_batch([position_counts], model)[0]


def detect_peaks_batch(
    histograms: List[Optional[Mapping[int, int]]], model=None
) -> List[List[Peak]]:
    """Detect peaks across multiple histograms with a single batched Cython/L1 GBDT predict call."""
    results: List[List[Peak]] = [[] for _ in histograms]
    summaries: Dict[int, RawHistogramSummary] = {}
    feature_rows: List[List[float]] = []
    row_owner: List[Tuple[int, dict]] = []

    for i, counts in enumerate(histograms):
        if not counts:
            continue
        hist = counts if isinstance(counts, TerminalHistogram) else TerminalHistogram(counts)
        summary = hist.summarize(compute_prominence=False)
        if summary is None:
            continue
        summaries[i] = summary
        if summary.peak_count == 0:
            pos = summary.start + summary.mode
            cnt = window_count(summary.padded, summary.mode, summary.mode)
            results[i].append(Peak(pos, cnt, pos, pos))
            continue

        for peak in rank_peaks_summary(summary):
            feature_rows.append(
                [
                    summary.var,
                    summary.skew,
                    float(summary.peak_count),
                    float(peak["peak_width"]),
                    summary.entropy,
                    summary.mean_height,
                    float(peak["peak_heights"]),
                    float(peak["relative_height"]),
                ]
            )
            row_owner.append((i, peak))

    if feature_rows:
        active_model = get_polya_model() if model is None else model
        feat_mat = np.asarray(feature_rows, dtype=np.float32)
        predicted = active_model.predict(feat_mat)
        for (i, peak), keep in zip(row_owner, predicted):
            if int(keep) != 1:
                continue
            summary = summaries[i]
            start = summary.start
            pos = start + int(peak["peak_location"])
            cnt = window_count(summary.padded, peak["peak_left"], peak["peak_right"])
            results[i].append(
                Peak(pos, cnt, start + int(peak["peak_left"]), start + int(peak["peak_right"]))
            )
    return results
