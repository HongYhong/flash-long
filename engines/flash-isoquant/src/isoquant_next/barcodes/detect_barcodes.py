"""Parallel barcode & UMI detection driver using pysam.FastxFile and WorkerPoolExecutor."""

from __future__ import annotations

import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pysam

from isoquant_next.barcodes.callers import (
    CurioBarcodeDetector,
    SharedMemoryStereoBarcodeDetector,
    SharedMemoryStereoSplittingBarcodeDetector,
    TenXBarcodeDetector,
    TenXSplittingBarcodeDetector,
    TenXv2BarcodeDetector,
    TenXv2SplittingBarcodeDetector,
    UniversalSingleMoleculeExtractor,
    VisiumHDBarcodeDetector,
)
from isoquant_next.barcodes.common import iter_fastx_chunks, load_barcodes
from isoquant_next.barcodes.molecule_structure import MoleculeStructure
from isoquant_next.barcodes.results import ExtractionReadStats, ExtractionResult, ReadStats
from isoquant_next.config.enums import IsoQuantMode
from isoquant_next.core.concurrency import WorkerPoolExecutor
from isoquant_next.core.error_codes import IsoQuantExitCode

logger = logging.getLogger("IsoQuant")

BARCODE_CALLING_MODES: Dict[IsoQuantMode, type] = {
    IsoQuantMode.tenX_v3: TenXBarcodeDetector,
    IsoQuantMode.tenX_v2: TenXv2BarcodeDetector,
    IsoQuantMode.tenX_v3_split: TenXSplittingBarcodeDetector,
    IsoQuantMode.tenX_v2_split: TenXv2SplittingBarcodeDetector,
    IsoQuantMode.curio: CurioBarcodeDetector,
    IsoQuantMode.stereoseq_nosplit: SharedMemoryStereoBarcodeDetector,
    IsoQuantMode.stereoseq: SharedMemoryStereoSplittingBarcodeDetector,
    IsoQuantMode.visium_5prime: TenXBarcodeDetector,
    IsoQuantMode.visium_hd: VisiumHDBarcodeDetector,
    IsoQuantMode.custom_sc: UniversalSingleMoleculeExtractor,
}

BARCODE_FILES_REQUIRED: Dict[IsoQuantMode, List[int]] = {
    IsoQuantMode.tenX_v3: [1],
    IsoQuantMode.tenX_v2: [1],
    IsoQuantMode.tenX_v3_split: [1],
    IsoQuantMode.tenX_v2_split: [1],
    IsoQuantMode.curio: [1, 2],
    IsoQuantMode.stereoseq_nosplit: [1],
    IsoQuantMode.stereoseq: [1],
    IsoQuantMode.visium_5prime: [1],
    IsoQuantMode.visium_hd: [2],
    IsoQuantMode.custom_sc: [0],
}


def stats_file_name(file_name: str) -> str:
    return f"{file_name}.stats"


def get_umi_length(isoquant_mode: IsoQuantMode) -> int:
    if isoquant_mode in (IsoQuantMode.tenX_v3, IsoQuantMode.tenX_v3_split, IsoQuantMode.visium_5prime):
        return 12
    if isoquant_mode in (
        IsoQuantMode.tenX_v2,
        IsoQuantMode.tenX_v2_split,
        IsoQuantMode.stereoseq,
        IsoQuantMode.stereoseq_nosplit,
    ):
        return 10
    if isoquant_mode == IsoQuantMode.curio:
        return 9
    if isoquant_mode == IsoQuantMode.visium_hd:
        return 12
    return 0


class SimpleReadStorage:
    """Bounded read sequence buffer."""

    def __init__(self) -> None:
        self.read_ids: List[str] = []
        self.sequences: List[str] = []

    def add(self, read_id: str, seq: str) -> None:
        self.read_ids.append(read_id)
        self.sequences.append(seq)

    def clear(self) -> None:
        self.read_ids.clear()
        self.sequences.clear()


class BarcodeCaller:
    """Streaming FASTQ/FASTA/BAM barcode detector with automatic context-managed resource cleanup."""

    def __init__(
        self,
        output_file_name: str,
        barcode_detector,
        header: bool = False,
        output_sequences: Optional[str] = None,
    ) -> None:
        self.output_file_name = output_file_name
        self.barcode_detector = barcode_detector
        self.output_file = open(output_file_name, "w", encoding="utf-8")
        if header:
            self.output_file.write(barcode_detector.header() + "\n")
        self.out_fasta = None
        if output_sequences:
            self.out_fasta = open(output_sequences, "w", encoding="utf-8")
        res_type = barcode_detector.result_type()
        self.read_stats = ExtractionReadStats() if res_type is ExtractionResult else ReadStats()

    def __enter__(self) -> "BarcodeCaller":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def _process_read_normal(self, read_id: str, read_sequence: str) -> None:
        res = self.barcode_detector.find_barcode_umi(read_id, read_sequence)
        self.read_stats.add_read(res)
        self.output_file.write(str(res) + "\n")
        if self.out_fasta is not None and res.is_valid():
            start = res.get_fasta_segment_start() if hasattr(res, "get_fasta_segment_start") else 0
            end = res.get_fasta_segment_end(len(read_sequence)) if hasattr(res, "get_fasta_segment_end") else len(read_sequence)
            seg = read_sequence[start:end]
            if seg:
                self.out_fasta.write(f">{read_id}\n{seg}\n")

    def _process_read_split(self, read_id: str, read_sequence: str) -> None:
        res = self.barcode_detector.find_barcode_umi(read_id, read_sequence)
        patterns = getattr(res, "detected_patterns", [res])
        if len(patterns) > 1:
            self.read_stats.add_custom_stats("Split reads", 1)
        for idx, pat in enumerate(patterns):
            self.read_stats.add_read(pat)
            if len(patterns) > 1:
                pat.read_id = f"{read_id}__part{idx}"
            self.output_file.write(str(pat) + "\n")
            if self.out_fasta is not None and pat.is_valid():
                start = pat.get_fasta_segment_start() if hasattr(pat, "get_fasta_segment_start") else 0
                end = pat.get_fasta_segment_end(len(read_sequence)) if hasattr(pat, "get_fasta_segment_end") else len(read_sequence)
                seg = read_sequence[start:end]
                if seg:
                    self.out_fasta.write(f">{pat.read_id}\n{seg}\n")

    def _process_read(self, read_id: str, read_sequence: str) -> None:
        if hasattr(self.barcode_detector, "detect_multiple"):
            self._process_read_split(read_id, read_sequence)
        else:
            res = self.barcode_detector.find_barcode_umi(read_id, read_sequence)
            if hasattr(res, "detected_patterns"):
                patterns = res.detected_patterns
                if len(patterns) > 1 and hasattr(self.read_stats, "add_custom_stats"):
                    self.read_stats.add_custom_stats("Split reads", 1)
                for idx, pat in enumerate(patterns):
                    self.read_stats.add_read(pat)
                    if len(patterns) > 1:
                        pat.read_id = f"{read_id}__{idx}"
                    self.output_file.write(str(pat) + "\n")
            else:
                self.read_stats.add_read(res)
                self.output_file.write(str(res) + "\n")

    def process_chunk(self, read_chunk: List[Tuple[str, str]]) -> None:
        for read_id, read_sequence in read_chunk:
            self._process_read(read_id, read_sequence)

    def _process_fastx(self, input_file: str) -> None:
        for batch in iter_fastx_chunks(input_file, chunk_size=1000):
            for read_id, seq, _qual in batch:
                self._process_read(read_id, seq)

    def _process_bam(self, input_file: str) -> None:
        with pysam.AlignmentFile(input_file, "rb", check_sq=False) as bam:
            for aln in bam.fetch(until_eof=True):
                if aln.is_secondary or aln.is_supplementary:
                    continue
                if aln.query_name and aln.query_sequence:
                    self._process_read(aln.query_name, aln.query_sequence)

    def process(self, input_file: str) -> None:
        low = input_file.lower()
        if low.endswith(".bam") or low.endswith(".sam"):
            self._process_bam(input_file)
        else:
            self._process_fastx(input_file)

    def dump_stats(self, file_name: Optional[str] = None) -> None:
        out_path = file_name or stats_file_name(self.output_file_name)
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(str(self.read_stats))

    def get_stats(self):
        return self.read_stats

    def close(self) -> None:
        if self.output_file is not None and not self.output_file.closed:
            self.output_file.close()
        if self.out_fasta is not None and not self.out_fasta.closed:
            self.out_fasta.close()


def fastx_file_chunk_reader(handler: str, chunk_size: int = 1000) -> Iterable[List[Tuple[str, str]]]:
    for batch in iter_fastx_chunks(handler, chunk_size=chunk_size):
        yield [(rid, seq) for rid, seq, _q in batch]


def bam_file_chunk_reader(handler: str, chunk_size: int = 1000) -> Iterable[List[Tuple[str, str]]]:
    chunk: List[Tuple[str, str]] = []
    with pysam.AlignmentFile(handler, "rb", check_sq=False) as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_secondary or aln.is_supplementary:
                continue
            if aln.query_name and aln.query_sequence:
                chunk.append((aln.query_name, aln.query_sequence))
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
    if chunk:
        yield chunk


def create_barcode_caller(args):
    mode = getattr(args, "mode", None)
    if mode not in BARCODE_CALLING_MODES:
        logger.critical("Unsupported barcode calling mode %s", mode)
        sys.exit(int(IsoQuantExitCode.INVALID_PARAMETER))

    detector_cls = BARCODE_CALLING_MODES[mode]
    if mode == IsoQuantMode.custom_sc:
        mdf_path = getattr(args, "molecule_description", None) or getattr(args, "molecule", None)
        if not mdf_path or not os.path.isfile(mdf_path):
            logger.critical("Missing molecule description file %s", mdf_path)
            sys.exit(int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND))
        with open(mdf_path, "r", encoding="utf-8") as fh:
            ms = MoleculeStructure(fh)
        return UniversalSingleMoleculeExtractor(ms)

    bc_files = list(getattr(args, "barcode_whitelist", None) or getattr(args, "barcodes", None) or [])
    req = BARCODE_FILES_REQUIRED[mode]
    if len(bc_files) not in req:
        logger.critical("Expected %s barcode files for mode %s, got %d", req, mode, len(bc_files))
        sys.exit(int(IsoQuantExitCode.INVALID_PARAMETER))

    if len(bc_files) == 1:
        barcodes = load_barcodes(bc_files[0])
        return detector_cls(barcodes)
    bc1 = load_barcodes(bc_files[0])
    bc2 = load_barcodes(bc_files[1])
    if mode == IsoQuantMode.visium_hd:
        return detector_cls([bc1, bc2])
    return detector_cls((bc1, bc2))


def process_chunk(
    barcode_detector,
    read_chunk: List[Tuple[str, str]],
    output_file: str,
    num: int,
    out_fasta: Optional[str] = None,
) -> Tuple[str, Optional[str], object]:
    chunk_tsv = f"{output_file}.{num}.tmp"
    chunk_fa = f"{out_fasta}.{num}.tmp" if out_fasta else None
    with BarcodeCaller(chunk_tsv, barcode_detector, header=False, output_sequences=chunk_fa) as caller:
        caller.process_chunk(read_chunk)
        stats = caller.get_stats()
    return chunk_tsv, chunk_fa, stats


def _concat_files(src_paths: List[str], dst_path: str, append: bool = False, skip_first_header_line: bool = False) -> None:
    mode = "ab" if append else "wb"
    with open(dst_path, mode) as dst:
        for i, src_p in enumerate(src_paths):
            if not os.path.exists(src_p):
                continue
            with open(src_p, "rb") as src:
                if skip_first_header_line and i > 0:
                    src.readline()
                shutil.copyfileobj(src, dst, 1024 * 1024)


def _process_single_file_in_parallel(
    input_file: str,
    output_tsv: str,
    out_fasta: Optional[str],
    args,
    barcode_detector,
) -> object:
    threads = max(1, int(getattr(args, "threads", 1)))
    reader = bam_file_chunk_reader(input_file) if input_file.lower().endswith((".bam", ".sam")) else fastx_file_chunk_reader(input_file)
    tmp_tsvs: List[str] = []
    tmp_fas: List[str] = []
    agg_stats = ExtractionReadStats() if barcode_detector.result_type() is ExtractionResult else ReadStats()

    with WorkerPoolExecutor(max_workers=threads) as pool:
        futures = [
            pool.submit(process_chunk, barcode_detector, chunk, output_tsv, idx, out_fasta)
            for idx, chunk in enumerate(reader)
        ]
        for fut in futures:
            tsv_part, fa_part, st = fut.result()
            tmp_tsvs.append(tsv_part)
            if fa_part:
                tmp_fas.append(fa_part)
            agg_stats.read_count += st.read_count
            agg_stats.bc_count += st.bc_count
            agg_stats.umi_count += st.umi_count
            if hasattr(agg_stats, "additional_attributes_counts") and hasattr(st, "additional_attributes_counts"):
                for k, v in st.additional_attributes_counts.items():
                    agg_stats.additional_attributes_counts[k] += v
            if hasattr(agg_stats, "pattern_counts") and hasattr(st, "pattern_counts"):
                for k, v in st.pattern_counts.items():
                    agg_stats.pattern_counts[k] += v

    with open(output_tsv, "w", encoding="utf-8") as dst:
        dst.write(barcode_detector.header() + "\n")
    _concat_files(tmp_tsvs, output_tsv, append=True)
    for p in tmp_tsvs:
        if os.path.exists(p):
            os.remove(p)
    if out_fasta and tmp_fas:
        if os.path.exists(out_fasta):
            os.remove(out_fasta)
        _concat_files(tmp_fas, out_fasta, append=False)
        for p in tmp_fas:
            if os.path.exists(p):
                os.remove(p)
    return agg_stats


def _resolve_barcode_io_lists(args) -> Tuple[List[str], List[str], List[Optional[str]]]:
    input_files = getattr(args, "input", None) or getattr(args, "fastq", None) or []
    input_list = [str(input_files)] if isinstance(input_files, (str, Path)) else [str(x) for x in input_files]

    raw_tsv = getattr(args, "output_tsv", None)
    if raw_tsv is None:
        out_tsvs = [str(args.output)]
    elif isinstance(raw_tsv, (str, Path)):
        out_tsvs = [str(raw_tsv)]
    else:
        out_tsvs = [str(x) for x in raw_tsv]

    raw_fa = getattr(args, "out_fasta", None)
    if raw_fa is None:
        out_fastas: List[Optional[str]] = [None] * len(input_list)
    elif isinstance(raw_fa, (str, Path)):
        out_fastas = [str(raw_fa)]
    else:
        out_fastas = [str(x) if x is not None else None for x in raw_fa]
    return input_list, out_tsvs, out_fastas


def process_single_thread(args) -> None:
    detector = create_barcode_caller(args)
    input_files, output_tsvs, out_fastas = _resolve_barcode_io_lists(args)
    try:
        if len(output_tsvs) == len(input_files):
            for i, inf in enumerate(input_files):
                out_fa = out_fastas[i] if i < len(out_fastas) else None
                with BarcodeCaller(output_tsvs[i], detector, header=True, output_sequences=out_fa) as caller:
                    caller.process(inf)
                    caller.dump_stats()
        else:
            out_fa = out_fastas[0] if out_fastas else None
            with BarcodeCaller(output_tsvs[0], detector, header=True, output_sequences=out_fa) as caller:
                for inf in input_files:
                    caller.process(inf)
                caller.dump_stats()
    finally:
        indexer = getattr(detector, "barcode_indexer", None)
        if indexer is not None and hasattr(indexer, "cleanup"):
            indexer.cleanup()


def process_in_parallel(args) -> None:
    detector = create_barcode_caller(args)
    input_files, output_tsvs, out_fastas = _resolve_barcode_io_lists(args)
    part_tsvs: List[str] = []
    part_fas: List[str] = []
    total_stats = None
    per_file_outputs = len(output_tsvs) == len(input_files)
    try:
        for i, inf in enumerate(input_files):
            p_tsv = output_tsvs[i] if per_file_outputs else (
                output_tsvs[0] if len(input_files) == 1 else f"{output_tsvs[0]}.sample{i}"
            )
            fa_target = out_fastas[i] if (per_file_outputs and i < len(out_fastas)) else (out_fastas[0] if out_fastas else None)
            p_fa = fa_target if (per_file_outputs or len(input_files) == 1 or not fa_target) else f"{fa_target}.sample{i}"
            st = _process_single_file_in_parallel(inf, p_tsv, p_fa, args, detector)
            part_tsvs.append(p_tsv)
            if p_fa:
                part_fas.append(p_fa)
            if per_file_outputs:
                with open(stats_file_name(p_tsv), "w", encoding="utf-8") as fh:
                    fh.write(str(st))
            else:
                total_stats = st if total_stats is None else st
        if not per_file_outputs and len(part_tsvs) > 1:
            _concat_files(part_tsvs, output_tsvs[0], append=False, skip_first_header_line=True)
            for p in part_tsvs:
                os.remove(p)
        if not per_file_outputs and out_fastas and out_fastas[0] and len(part_fas) > 1:
            _concat_files(part_fas, out_fastas[0], append=False)
            for p in part_fas:
                os.remove(p)
        if not per_file_outputs and total_stats is not None:
            with open(stats_file_name(output_tsvs[0]), "w", encoding="utf-8") as fh:
                fh.write(str(total_stats))
    finally:
        indexer = getattr(detector, "barcode_indexer", None)
        if indexer is not None and hasattr(indexer, "cleanup"):
            indexer.cleanup()
