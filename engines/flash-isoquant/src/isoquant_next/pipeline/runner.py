"""End-to-end Stage Graph pipeline runner wiring barcodes, `.iqidx`, minimap2, discovery/quantification, and fusion."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
from typing import Dict, List, Optional

import pysam

from isoquant_next.annotation.gtf_validator import check_gtf_duplicates_streaming
from isoquant_next.annotation.iqidx import (
    compile_gffutils_db_to_iqidx,
    compile_gtf_to_iqidx,
)
from isoquant_next.barcodes.detect_barcodes import process_in_parallel, process_single_thread
from isoquant_next.core.concurrency import WorkerPoolExecutor
from isoquant_next.core.error_codes import IsoQuantError, IsoQuantExitCode
from isoquant_next.fusion import FusionDetector
from isoquant_next.mapping.read_mapper import ReadMapper, index_bam_safe
from isoquant_next.pipeline.dataset_processor import DatasetProcessor
from isoquant_next.pipeline.input_data import InputDataType, SampleData
from isoquant_next.pipeline.runtime_params import iqidx_to_bed12
from isoquant_next.pipeline.stage_graph import StageGraphExecutor, StageNode, StageStateTracker
from isoquant_next.quantification import combine_counts

logger = logging.getLogger("IsoQuant")


@dataclass
class _BarcodeCallingArgs:
    input: List[str]
    barcodes: Optional[List[str]]
    mode: object
    output_tsv: List[str]
    out_fasta: Optional[List[str]]
    tmp_dir: str
    threads: int
    molecule: Optional[str] = None


def _run_barcode_worker(bc_args: _BarcodeCallingArgs) -> None:
    if bc_args.threads == 1:
        process_single_thread(bc_args)
    else:
        process_in_parallel(bc_args)


def _prepare_reference_genome(args) -> None:
    if not getattr(args, "reference", None):
        return
    ref_path = args.reference
    ref_dir = os.path.dirname(ref_path) or "."
    ref_base = os.path.basename(ref_path)
    stem, ext = os.path.splitext(ref_base)

    args.fai_file_name = f"{ref_path}.fai"
    if not os.path.exists(args.fai_file_name) and not os.access(ref_dir, os.W_OK):
        args.fai_file_name = os.path.join(args.output, f"{ref_base}.fai")

    if ext.lower() in (".gz", ".gzip", ".bgz"):
        gunzipped = os.path.join(args.output, stem)
        if not os.path.exists(gunzipped) or not args.resume:
            with gzip.open(ref_path, "rt", encoding="utf-8", errors="replace") as src, open(gunzipped, "w", encoding="utf-8") as dst:
                shutil.copyfileobj(src, dst)
        args.reference = gunzipped
        args.fai_file_name = f"{gunzipped}.fai"


def _prepare_iqidx_annotation(args) -> None:
    if not getattr(args, "genedb", None):
        return
    src = Path(args.genedb)
    args.original_annotation = str(src)
    if str(src).lower().endswith(".iqidx"):
        return

    out_dir = Path(getattr(args, "genedb_output", None) or args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    target_iqidx = out_dir / f"{src.stem}.iqidx"
    if (
        not args.clean_start
        and target_iqidx.exists()
        and target_iqidx.stat().st_mtime >= src.stat().st_mtime
    ):
        args.genedb = str(target_iqidx)
        return

    if str(src).lower().endswith(".db"):
        compile_gffutils_db_to_iqidx(src, target_iqidx)
    else:
        gtf_to_compile = src
        tmp_clean: Optional[str] = None
        if args.gtf_check:
            clean_out = out_dir / f"{src.stem}.cleaned.gtf"
            was_mod, clean_path = check_gtf_duplicates_streaming(src, clean_out)
            if was_mod and clean_path:
                gtf_to_compile = Path(clean_path)
                tmp_clean = clean_path
            elif clean_out.exists():
                clean_out.unlink()
        try:
            compile_gtf_to_iqidx(gtf_to_compile, target_iqidx)
        finally:
            if tmp_clean and os.path.exists(tmp_clean):
                os.remove(tmp_clean)
    args.genedb = str(target_iqidx)


def _call_barcodes_stage(args) -> None:
    if getattr(args, "barcoded_bam", False):
        return
    if getattr(args, "barcoded_reads", None):
        args.input_data.samples[0].barcoded_reads = list(args.barcoded_reads)
        return

    for sample in args.input_data.samples:
        os.makedirs(sample.out_dir, exist_ok=True)
        os.makedirs(sample.aux_dir, exist_ok=True)
        input_files = [files[0] for files in sample.file_list]
        out_tsvs = [f"{sample.barcodes_tsv}_{i}.tsv" for i in range(len(input_files))]
        done_markers = [f"{sample.barcodes_done}_{i}.tsv" for i in range(len(input_files))]
        out_fastas = (
            [f"{sample.split_reads_fasta}_{i}.fa" for i in range(len(input_files))]
            if args.mode.produces_new_fasta()
            else None
        )
        if args.resume and all(os.path.exists(d) for d in done_markers):
            sample.barcoded_reads.extend(out_tsvs)
            if out_fastas is not None:
                sample.file_list = [[fa] for fa in out_fastas]
            continue

        for d in done_markers:
            if os.path.exists(d):
                os.remove(d)
        bc_threads = 1 if args.mode.enforces_single_thread() else args.threads
        bc_args = _BarcodeCallingArgs(
            input_files, args.barcode_whitelist, args.mode, out_tsvs, out_fastas, sample.aux_dir, bc_threads,
            molecule=getattr(args, "molecule", None),
        )
        with WorkerPoolExecutor(max_workers=1, debug=getattr(args, "debug", False), log_file=getattr(args, "log_file", None)) as pool:
            pool.submit(_run_barcode_worker, bc_args).result()

        for d in done_markers:
            open(d, "w", encoding="utf-8").close()
        sample.barcoded_reads.extend(out_tsvs)
        if out_fastas is not None:
            sample.file_list = [[fa] for fa in out_fastas]


def _load_align_cache(path: str) -> dict:
    if not os.path.isfile(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}


def _save_align_cache(path: str, payload: dict) -> None:
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    os.replace(tmp, path)


def _map_reads_stage(args) -> None:
    if not args.input_data.input_type.needs_mapping():
        return

    mapper = ReadMapper(args.data_type, threads=args.threads)
    junc_bed = None
    if not args.no_junc_bed and args.genedb:
        if getattr(args, "junc_bed_file", None) and os.path.exists(args.junc_bed_file):
            junc_bed = os.path.abspath(args.junc_bed_file)
        else:
            bed_path = os.path.join(args.output, f"{Path(args.genedb).stem}.junc.bed")
            if args.clean_start or not os.path.exists(bed_path):
                iqidx_to_bed12(args.genedb, bed_path)
            junc_bed = bed_path

    align_target = getattr(args, "index", None) or args.reference
    stranded = str(getattr(args, "stranded", "none") or "none")
    mapping_opts = getattr(args, "mapping_options", None)
    cache = {} if args.clean_start else _load_align_cache(args.alignment_config_path)
    new_samples: List[SampleData] = []
    for sample in args.input_data.samples:
        os.makedirs(sample.aux_dir, exist_ok=True)
        readable: Dict[str, str] = {}
        bam_files: List[List[str]] = []
        for lib in sample.file_list:
            for reads_file in lib:
                abs_reads = os.path.abspath(reads_file)
                cache_key = f"{abs_reads}::{os.path.abspath(align_target)}::{args.data_type}::{junc_bed or ''}"
                cached_entry = cache.get(cache_key)
                bam_out: Optional[str] = None
                if (
                    cached_entry
                    and os.path.exists(cached_entry.get("bam", ""))
                    and os.path.getmtime(abs_reads) == cached_entry.get("reads_mtime")
                ):
                    bam_out = cached_entry["bam"]
                else:
                    if args.input_data.input_type == InputDataType.unmapped_bam:
                        fq_tmp = os.path.join(sample.aux_dir, f"{Path(abs_reads).stem}.unmapped.fastq")
                        if not mapper.samtools_bin:
                            raise IsoQuantError(IsoQuantExitCode.SAMTOOLS_FAILED, "samtools executable not found")
                        subprocess.run([mapper.samtools_bin, "fastq", abs_reads], check=True, stdout=open(fq_tmp, "wb"))
                        mapped_p = mapper.align_to_sorted_bam(
                            fq_tmp, align_target, sample.aux_dir, junc_bed=junc_bed, sample_label=sample.prefix,
                            stranded=stranded, mapping_options=mapping_opts,
                        )
                        if os.path.exists(fq_tmp):
                            os.remove(fq_tmp)
                    else:
                        mapped_p = mapper.align_to_sorted_bam(
                            abs_reads, align_target, sample.aux_dir, junc_bed=junc_bed, sample_label=sample.prefix,
                            stranded=stranded, mapping_options=mapping_opts,
                        )
                    bam_out = str(mapped_p)
                    index_bam_safe(bam_out, threads=args.threads)
                    cache[cache_key] = {"bam": bam_out, "reads_mtime": os.path.getmtime(abs_reads)}
                    _save_align_cache(args.alignment_config_path, cache)

                bam_files.append([bam_out])
                if reads_file in sample.readable_names_dict:
                    readable[bam_out] = sample.readable_names_dict[reads_file]

        new_samples.append(
            SampleData(bam_files, sample.prefix, sample.out_dir, readable, sample.illumina_bam, sample.barcoded_reads)
        )

    args.input_data.samples = new_samples
    args.input_data.input_type = InputDataType.bam


def _detect_fusions_stage(args) -> None:
    bam_files = [f for s in args.input_data.samples for lib in s.file_list for f in lib if os.path.isfile(f)]
    if not args.genedb or not bam_files:
        return
    with FusionDetector(bam_files[0], args.genedb, reference_fasta=args.reference) as fd:
        for sample in args.input_data.samples:
            sample_bams = [f for lib in sample.file_list for f in lib if os.path.isfile(f)]
            if not sample_bams:
                continue
            fd.clear_state()
            for bp in sample_bams:
                fd.bam_path = bp
                fd.detect_fusions()
            fd.report(output_path=os.path.join(sample.out_dir, f"{sample.prefix}.fusions.tsv"))


def run_pipeline(args, *, config_fingerprint: str = "default") -> List[str]:
    """Execute full IsoQuant Stage Graph DAG deterministically with atomic state tracking."""
    os.makedirs(args.output, exist_ok=True)
    tracker = StageStateTracker(os.path.join(args.output, ".stage_state.json"), fingerprint=config_fingerprint)

    stages = [
        StageNode(
            "barcode_calling",
            lambda: _call_barcodes_stage(args),
            enabled=args.mode.needs_barcode_calling(),
            on_skip=lambda: _call_barcodes_stage(args),
        ),
        StageNode(
            "reference_and_annotation",
            lambda: (_prepare_reference_genome(args), _prepare_iqidx_annotation(args)),
            enabled=True,
            on_skip=lambda: (_prepare_reference_genome(args), _prepare_iqidx_annotation(args)),
        ),
        StageNode(
            "read_mapping",
            lambda: _map_reads_stage(args),
            deps=("barcode_calling", "reference_and_annotation"),
            enabled=args.input_data.input_type.needs_mapping(),
            on_skip=lambda: _map_reads_stage(args),
        ),
        StageNode(
            "dataset_processing",
            lambda: DatasetProcessor(args).process_all_samples(args.input_data),
            deps=("read_mapping", "reference_and_annotation"),
            enabled=not getattr(args, "run_aligner_only", False),
        ),
        StageNode(
            "combine_counts",
            lambda: combine_counts(args.input_data, args.output),
            deps=("dataset_processing",),
            enabled=(
                not getattr(args, "run_aligner_only", False)
                and len(args.input_data.samples) > 1
                and bool(args.genedb)
                and args.run_quantification
            ),
        ),
        StageNode(
            "fusion_detection",
            lambda: _detect_fusions_stage(args),
            deps=("read_mapping", "reference_and_annotation", "dataset_processing"),
            enabled=not getattr(args, "run_aligner_only", False) and bool(getattr(args, "fusion", False)),
        ),
    ]
    executor = StageGraphExecutor(stages, tracker=tracker)
    executed = executor.execute(resume=bool(getattr(args, "resume", False)))
    logger.info("IsoQuant pipeline completed stages: %s", ", ".join(executed) if executed else "all cached")
    return executed
