"""Terminal (polyA / TSS) prediction counters using sparse TerminalHistogram, L1TreeEnsemble, and np.searchsorted."""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from isoquant_next.assignment.types import ReadAssignment, ReadAssignmentType
from isoquant_next.terminal.peaks import (
    ANNOTATION_TOLERANCE,
    FEATURE_COLUMNS,
    HISTOGRAM_PAD,
    POLYA_MODEL_PATH,
    TSS_MODEL_PATH,
    RawHistogramSummary,
    TerminalHistogram,
    rank_peaks_summary,
    window_count,
)
from isoquant_next.terminal.tree_model import L1TreeEnsemble

logger = logging.getLogger("IsoQuant")

_ACCEPTED_ASSIGNMENT_TYPES = frozenset(
    {
        ReadAssignmentType.unique,
        ReadAssignmentType.unique_minor_difference,
        ReadAssignmentType.inconsistent,
        ReadAssignmentType.inconsistent_non_intronic,
    }
)

EMPTY_COLUMNS = ["chromosome", "transcript_id", "gene_id", "prediction", "counts", "flag"]
EMPTY_COLUMNS_GROUPED = EMPTY_COLUMNS + ["counts_byGroup", "group_id"]
TRAINING_COLUMNS = FEATURE_COLUMNS + ["chromosome", "true_peak"]
TRAINING_SUFFIX = ".training.csv"


class _TranscriptEntry:
    __slots__ = ("chr", "gene_id", "annotated", "histogram", "group_positions")

    def __init__(self, chrom: str, gene_id: str, annotated: int) -> None:
        self.chr = chrom
        self.gene_id = gene_id
        self.annotated = annotated
        self.histogram = TerminalHistogram()
        self.group_positions: Dict[int, List[int]] = defaultdict(list)


class TerminalCounter:
    """Zero-pandas, zero-xgboost polyA/TSS terminal peak counter with np.searchsorted grouped counting."""

    TRAINING_ARG: str = ""

    def __init__(
        self,
        args,
        output_prefix: str,
        model_path: Union[str, Path],
        string_pools=None,
        group_index: int = 0,
    ) -> None:
        self.ignore_read_groups = string_pools is None
        self.output_prefix = output_prefix
        self.output_file = output_prefix
        self.output_counts_file_name = output_prefix
        self.output_tpm_file_name = None
        self.output_stats_file_name = None
        self.usable_file_name = None
        open(self.output_file, "w").close()

        self.args = args
        self.string_pools = string_pools
        self.group_index = group_index
        self._model_path = Path(model_path)
        self._model = None

        self._collecting_training = bool(
            self.TRAINING_ARG and getattr(args, self.TRAINING_ARG, None)
        )
        self._training_csv_path = (
            output_prefix + TRAINING_SUFFIX if self._collecting_training else None
        )
        if self._training_csv_path:
            open(self._training_csv_path, "w").close()

        self.transcripts: Dict[str, _TranscriptEntry] = {}
        self._all_transcripts: Dict[str, _TranscriptEntry] = {}
        self.last_gene_predictions: List[int] = []

    @property
    def model(self) -> L1TreeEnsemble:
        if self._model is None:
            self._model = L1TreeEnsemble(self._model_path)
        return self._model

    def add_confirmed_features(self, features) -> None:
        return

    def add_unassigned(self, read_assignment) -> None:
        return

    def add_unaligned(self, n_reads: int = 1) -> None:
        return

    def finalize(self, args=None) -> None:
        return

    def get_output_file_handler(self):
        return open(self.output_file, "a")

    def _passes_filter(self, read_assignment: ReadAssignment) -> bool:
        raise NotImplementedError

    def _extract_position(self, read_assignment: ReadAssignment) -> Optional[int]:
        raise NotImplementedError

    def _annotated_position(self, read_assignment: ReadAssignment, transcript_id: str) -> int:
        raise NotImplementedError

    def _read_group_id(self, read_assignment: ReadAssignment) -> int:
        if not read_assignment.read_group_ids:
            return 0
        return read_assignment.read_group_ids[self.group_index]

    def add_read_info(self, read_assignment: Optional[ReadAssignment] = None) -> None:
        if read_assignment is None:
            return
        if read_assignment.assignment_type not in _ACCEPTED_ASSIGNMENT_TYPES:
            return
        if not self._passes_filter(read_assignment):
            return
        if not read_assignment.corrected_exons or not read_assignment.isoform_matches:
            return
        position = self._extract_position(read_assignment)
        if position is None:
            return

        isoform_match = read_assignment.isoform_matches[0]
        transcript_id = isoform_match.assigned_transcript
        if transcript_id is None:
            return

        entry = self.transcripts.get(transcript_id)
        if entry is None:
            entry = _TranscriptEntry(
                read_assignment.chr_id,
                isoform_match.assigned_gene,
                self._annotated_position(read_assignment, transcript_id),
            )
            self.transcripts[transcript_id] = entry

        pos = int(position)
        entry.histogram.add(pos, 1)
        if not self.ignore_read_groups:
            gid = self._read_group_id(read_assignment)
            entry.group_positions[gid].append(pos)

    def flush(self) -> None:
        if not self.ignore_read_groups or self._collecting_training:
            return
        rows = self._predict_rows()
        self.last_gene_predictions = [r["prediction"] for r in rows] if rows else []
        self._merge_into_all()
        self.transcripts = {}

    def _merge_into_all(self) -> None:
        for tid, entry in self.transcripts.items():
            existing = self._all_transcripts.get(tid)
            if existing is None:
                self._all_transcripts[tid] = entry
            else:
                existing.histogram.merge(entry.histogram)

    def dump(self) -> None:
        if self._collecting_training:
            self._dump_training()
            return
        if not self.ignore_read_groups:
            rows = self._predict_rows()
            if not rows:
                self._write_empty()
            else:
                self._write_grouped(rows)
            return

        self.flush()
        self.transcripts = self._all_transcripts
        rows = self._predict_rows()
        if not rows:
            self._write_empty()
        else:
            self._write_ungrouped(rows)

    def _predict_rows(self) -> List[dict]:
        if not self.transcripts:
            return []

        zero_rows: List[dict] = []
        cand_features: List[List[float]] = []
        cand_meta: List[Tuple[str, _TranscriptEntry, RawHistogramSummary, dict]] = []

        for tid, entry in self.transcripts.items():
            summary = entry.histogram.summarize(compute_prominence=False)
            if summary is None:
                continue
            if summary.peak_count == 0:
                pred = summary.start + summary.mode
                cnt = window_count(summary.padded, summary.mode, summary.mode)
                flag = "Novel" if abs(pred - entry.annotated) > ANNOTATION_TOLERANCE else "Known"
                zero_rows.append(
                    {
                        "chromosome": entry.chr,
                        "transcript_id": tid,
                        "gene_id": entry.gene_id,
                        "prediction": pred,
                        "counts": cnt,
                        "flag": flag,
                        "low_abs": pred,
                        "high_abs": pred,
                    }
                )
                continue

            for peak in rank_peaks_summary(summary):
                cand_features.append(
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
                cand_meta.append((tid, entry, summary, peak))

        kept_rows: List[dict] = []
        if cand_features:
            predicted = self.model.predict(np.asarray(cand_features, dtype=np.float32))
            for (tid, entry, summary, peak), keep in zip(cand_meta, predicted):
                if int(keep) != 1:
                    continue
                pred = summary.start + int(peak["peak_location"])
                cnt = window_count(summary.padded, peak["peak_left"], peak["peak_right"])
                flag = "Novel" if abs(pred - entry.annotated) > ANNOTATION_TOLERANCE else "Known"
                kept_rows.append(
                    {
                        "chromosome": entry.chr,
                        "transcript_id": tid,
                        "gene_id": entry.gene_id,
                        "prediction": pred,
                        "counts": cnt,
                        "flag": flag,
                        "low_abs": summary.start + int(peak["peak_left"]),
                        "high_abs": summary.start + int(peak["peak_right"]),
                    }
                )

        return zero_rows + kept_rows

    def _dump_training(self) -> None:
        if not self.transcripts:
            self._write_empty()
            return
        with open(self._training_csv_path, "w", encoding="utf-8") as fh:
            fh.write(",".join(TRAINING_COLUMNS) + "\n")
            zero_lines: List[str] = []
            peak_lines: List[str] = []
            for tid, entry in self.transcripts.items():
                summary = entry.histogram.summarize(compute_prominence=True)
                if summary is None:
                    continue
                if summary.peak_count == 0:
                    pred = summary.start + summary.mode
                    true_pk = 1 if abs(pred - entry.annotated) <= ANNOTATION_TOLERANCE else 0
                    h_mode = int(summary.padded[summary.mode + HISTOGRAM_PAD])
                    vals = [
                        summary.var,
                        summary.skew,
                        0,
                        "",
                        summary.entropy,
                        summary.mean_height,
                        h_mode,
                        "",
                        entry.chr,
                        true_pk,
                    ]
                    zero_lines.append(",".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in vals))
                else:
                    for peak in rank_peaks_summary(summary):
                        pred = summary.start + int(peak["peak_location"])
                        true_pk = 1 if abs(pred - entry.annotated) <= ANNOTATION_TOLERANCE else 0
                        vals = [
                            summary.var,
                            summary.skew,
                            summary.peak_count,
                            peak["peak_width"],
                            summary.entropy,
                            summary.mean_height,
                            peak["peak_heights"],
                            peak["relative_height"],
                            entry.chr,
                            true_pk,
                        ]
                        peak_lines.append(",".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in vals))
            for line in zero_lines + peak_lines:
                fh.write(line + "\n")
        self._write_empty()

    def _write_empty(self) -> None:
        cols = EMPTY_COLUMNS if self.ignore_read_groups else EMPTY_COLUMNS_GROUPED
        with open(self.output_file, "w", encoding="utf-8") as fh:
            fh.write("\t".join(cols) + "\n")

    def _write_ungrouped(self, rows: List[dict]) -> None:
        with open(self.output_file, "w", encoding="utf-8") as fh:
            fh.write("\t".join(EMPTY_COLUMNS) + "\n")
            for r in rows:
                fh.write(
                    f"{r['chromosome']}\t{r['transcript_id']}\t{r['gene_id']}\t"
                    f"{r['prediction']}\t{r['counts']}\t{r['flag']}\n"
                )

    def _write_grouped(self, rows: List[dict]) -> None:
        seen_groups = sorted(
            {gid for entry in self.transcripts.values() for gid in entry.group_positions.keys()}
        )
        if not seen_groups:
            self._write_empty()
            return

        # Pre-sort group arrays per (transcript_id, group_id) once for O(log N) np.searchsorted window queries.
        sorted_cache: Dict[str, Dict[int, np.ndarray]] = {}
        for tid, entry in self.transcripts.items():
            gmap: Dict[int, np.ndarray] = {}
            for gid, pos_list in entry.group_positions.items():
                arr = np.asarray(pos_list, dtype=np.int64)
                arr.sort()
                gmap[gid] = arr
            sorted_cache[tid] = gmap

        group_names = {gid: self._group_name(gid) for gid in seen_groups}
        with open(self.output_file, "w", encoding="utf-8") as fh:
            fh.write("\t".join(EMPTY_COLUMNS_GROUPED) + "\n")
            for r in rows:
                gmap = sorted_cache[r["transcript_id"]]
                low = r["low_abs"]
                high = r["high_abs"]
                for gid in seen_groups:
                    arr = gmap.get(gid)
                    if arr is None or arr.size == 0:
                        continue
                    count = int(np.searchsorted(arr, high, side="right") - np.searchsorted(arr, low, side="left"))
                    if count == 0:
                        continue
                    fh.write(
                        f"{r['chromosome']}\t{r['transcript_id']}\t{r['gene_id']}\t"
                        f"{r['prediction']}\t{r['counts']}\t{r['flag']}\t{count}\t{group_names[gid]}\n"
                    )

    def _group_name(self, group_id: int) -> str:
        if self.string_pools is None:
            return "NA"
        pool = self.string_pools.get_read_group_pool(self.group_index)
        return pool.get_str(group_id)


class PolyACounter(TerminalCounter):
    """Predicts polyA cleavage sites per transcript."""

    TRAINING_ARG = "collect_polya_training"

    def __init__(self, args, output_prefix: str, string_pools=None, group_index: int = 0) -> None:
        super().__init__(args, output_prefix, POLYA_MODEL_PATH, string_pools=string_pools, group_index=group_index)

    def _passes_filter(self, read_assignment: ReadAssignment) -> bool:
        return bool(read_assignment.polyA_found and read_assignment.polya_info)

    def _extract_position(self, read_assignment: ReadAssignment) -> Optional[int]:
        info = read_assignment.polya_info
        if read_assignment.strand == "+" and info.external_polya_pos != -1:
            return info.external_polya_pos
        if read_assignment.strand == "-" and info.external_polyt_pos != -1:
            return info.external_polyt_pos
        return None

    def _annotated_position(self, read_assignment: ReadAssignment, transcript_id: str) -> int:
        exons = read_assignment.gene_info.all_isoforms_exons[transcript_id]
        if read_assignment.strand == "+":
            return exons[-1][1] + 1
        return exons[0][0] - 1


class TSSCounter(TerminalCounter):
    """Predicts transcription start sites per transcript."""

    TRAINING_ARG = "collect_tss_training"

    def __init__(self, args, output_prefix: str, string_pools=None, group_index: int = 0) -> None:
        super().__init__(args, output_prefix, TSS_MODEL_PATH, string_pools=string_pools, group_index=group_index)

    def _passes_filter(self, read_assignment: ReadAssignment) -> bool:
        return read_assignment.strand in ("+", "-")

    def _extract_position(self, read_assignment: ReadAssignment) -> Optional[int]:
        exons = read_assignment.corrected_exons
        if not exons:
            return None
        if read_assignment.strand == "+":
            return exons[0][0]
        return exons[-1][1]

    def _annotated_position(self, read_assignment: ReadAssignment, transcript_id: str) -> int:
        exons = read_assignment.gene_info.all_isoforms_exons[transcript_id]
        if read_assignment.strand == "+":
            return exons[0][0] - 1
        return exons[-1][1] + 1
