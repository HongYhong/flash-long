"""Main CLI entrypoint `isoquant-next` with JSON manifest persistence and atomic `--resume`."""

from __future__ import annotations

import argparse
import glob
import json
import logging
import os
from pathlib import Path
import shutil
import sys
from typing import Any, Dict, List, Optional, Sequence

import pysam

from dataclasses import replace

from isoquant_next import __version__
from isoquant_next.alignment.read_groups import parse_barcode2spot_spec
from isoquant_next.config.enums import (
    ANALYSIS_ALIASES,
    ANALYSIS_CHOICES,
    AnalysisType,
    COUNTING_STRATEGIES,
    CountingStrategy,
    DATA_TYPE_ALIASES,
    GroupedOutputFormat,
    ISOQUANT_MODES,
    IsoQuantMode,
    NormalizationMethod,
    PolyATrimmed,
    PolyAUsageStrategies,
    StrandnessReportingLevel,
)
from isoquant_next.config.settings import PipelineConfig, QuantificationConfig
from isoquant_next.core.error_codes import IsoQuantError, IsoQuantExitCode
from isoquant_next.core.logging_setup import configure_logging, pin_single_thread_math_env
from isoquant_next.core.utils import LARGE_OUTPUT_TYPES
from isoquant_next.pipeline.input_data import InputDataStorage, InputDataType
from isoquant_next.pipeline.runner import run_pipeline
from isoquant_next.pipeline.runtime_params import PipelineParamsView

logger = logging.getLogger("IsoQuant")

SUPPORTED_ALIGNERS = ("starlong", "minimap2")
SUPPORTED_STRANDEDNESS = ("forward", "reverse", "none")

_CLI_OPTIONS_FILE = "run_cli_options.json"
_MANIFEST_FILE = "run_manifest.json"


def _bool_str(s: str) -> bool:
    low = s.lower()
    if low in ("true", "yes", "1", "y", "t"):
        return True
    if low in ("false", "no", "0", "n", "f"):
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got {s!r}")


def _abspath_opt(val: Optional[str]) -> Optional[str]:
    return os.path.abspath(val) if val else None


def _abspath_list(vals: Optional[Sequence[str]]) -> Optional[List[str]]:
    return [os.path.abspath(v) for v in vals] if vals is not None else None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="isoquant-next", formatter_class=argparse.RawDescriptionHelpFormatter)
    ref_g = parser.add_argument_group("Reference data")
    in_g = parser.add_argument_group("Input data")
    out_g = parser.add_argument_group("Output naming")
    pipe_g = parser.add_argument_group("Pipeline options")
    algo_g = parser.add_argument_group("Algorithm settings")
    setup_g = parser.add_argument_group("Output configuration")
    aln_g = parser.add_argument_group("Aligner settings")
    flt_g = parser.add_argument_group("Read filtering options")
    sc_g = parser.add_argument_group("Single-cell/spatial-related options")

    parser.add_argument("--full_help", action="help", help="show full list of options")
    parser.add_argument("--test", action="store_true", default=False, help="run IsoQuant on toy dataset")
    parser.add_argument("--debug", action="store_true", default=False, help="debug log output")
    parser.add_argument("--version", "-v", action="version", version=f"IsoQuant {__version__}")

    ref_g.add_argument("--reference", "-r", type=str, help="reference genome in FASTA format")
    ref_g.add_argument("--genedb", "-g", type=str, help="gene annotation in GTF/GFF3/.db/.iqidx format")
    ref_g.add_argument("--complete_genedb", action="store_true", default=False)
    ref_g.add_argument("--discard_chr", nargs="+", type=str, default=[])
    ref_g.add_argument("--process_only_chr", nargs="+", type=str, default=None)
    ref_g.add_argument("--index", type=str, default=None)

    in_g.add_argument("--data_type", "-d", type=str, choices=list(DATA_TYPE_ALIASES.keys()))
    src_g = in_g.add_mutually_exclusive_group()
    src_g.add_argument("--bam", nargs="+", type=str)
    src_g.add_argument("--fastq", nargs="+", type=str)
    src_g.add_argument("--unmapped_bam", nargs="+", type=str)
    src_g.add_argument("--yaml", type=str)

    in_g.add_argument("--illumina_bam", nargs="+", type=str, default=None)
    in_g.add_argument("--read_group", nargs="+", type=str, default=None)
    in_g.add_argument("--read_assignments", nargs="+", type=str, default=None)
    in_g.add_argument("--stranded", type=str, choices=SUPPORTED_STRANDEDNESS, default="none")
    in_g.add_argument(
        "--polya_trimmed", type=str, default=PolyATrimmed.none.name, choices=[e.name for e in PolyATrimmed]
    )
    in_g.add_argument("--fl_data", action="store_true", default=False)

    out_g.add_argument("--output", "-o", type=str, default="isoquant_output")
    out_g.add_argument("--prefix", "-p", type=str, default="OUT")
    out_g.add_argument("--labels", "-l", nargs="+", type=str, default=None)

    pipe_g.add_argument("--analysis", nargs="+", type=str, choices=ANALYSIS_CHOICES, default=None)
    pipe_g.add_argument("--threads", "-t", type=int, default=16)
    pipe_g.add_argument("--count_exons", action="store_true", default=False)
    pipe_g.add_argument("--count_intron_retentions", action="store_true", default=False)
    pipe_g.add_argument("--no_model_construction", action="store_true", default=False)
    res_g = pipe_g.add_mutually_exclusive_group()
    res_g.add_argument("--resume", action="store_true", default=False)
    res_g.add_argument("--force", action="store_true", default=False)
    pipe_g.add_argument("--clean_start", action="store_true", default=False)
    pipe_g.add_argument("--run_aligner_only", action="store_true", default=False)
    pipe_g.add_argument("--no_gtf_check", dest="gtf_check", action="store_false", default=True)
    pipe_g.add_argument("--high_memory", action="store_true", default=False)
    pipe_g.add_argument("--keep_tmp", action="store_true", default=False)

    setup_g.add_argument("--check_canonical", action="store_true", default=False)
    setup_g.add_argument("--sqanti_output", action="store_true", default=False)
    setup_g.add_argument("--bam_tags", type=str, default=None)
    setup_g.add_argument(
        "--large_output",
        nargs="+",
        type=str,
        choices=list(LARGE_OUTPUT_TYPES),
        default=["read_info", "read2transcripts"],
    )
    setup_g.add_argument("--genedb_output", type=str, default=None)
    setup_g.add_argument("--no_gzip", dest="gzipped", action="store_false", default=True)
    setup_g.add_argument(
        "--normalization_method",
        type=str,
        choices=[e.name for e in NormalizationMethod],
        default=NormalizationMethod.simple.name,
    )
    setup_g.add_argument(
        "--counts_format",
        type=str,
        nargs="+",
        choices=[e.name for e in GroupedOutputFormat],
        default=[GroupedOutputFormat.default.name],
    )
    setup_g.add_argument("--emit_read_ids", action="store_true", default=False)
    setup_g.add_argument("--old_exon_count_format", action="store_true", default=False)

    aln_g.add_argument("--aligner", choices=SUPPORTED_ALIGNERS, type=str, default=None)
    aln_g.add_argument("--no_junc_bed", action="store_true", default=False)
    aln_g.add_argument("--junc_bed_file", type=str, default=None)
    aln_g.add_argument("--indexing_options", type=str, default=None)
    aln_g.add_argument("--mapping_options", type=str, default=None)

    flt_g.add_argument("--use_secondary", action="store_true", default=False)
    flt_g.add_argument("--min_mapq", type=int, default=None)
    flt_g.add_argument("--inconsistent_mapq_cutoff", type=int, default=5)
    flt_g.add_argument("--simple_alignments_mapq_cutoff", type=int, default=1)
    flt_g.add_argument("--max_coverage_small_chr", type=int, default=1_000_000)
    flt_g.add_argument("--max_coverage_normal_chr", type=int, default=-1)

    sc_g.add_argument("--mode", "-m", type=str, choices=ISOQUANT_MODES, default=IsoQuantMode.bulk.name)
    sc_g.add_argument("--barcode_whitelist", type=str, nargs="+", default=None)
    sc_g.add_argument("--barcoded_reads", type=str, nargs="+", default=None)
    sc_g.add_argument("--barcoded_bam", action="store_true", default=False)
    sc_g.add_argument("--barcode2spot", type=str, default=None)
    sc_g.add_argument("--molecule", type=str, default=None)
    sc_g.add_argument("--barcode_tag", type=str, default="CB")
    sc_g.add_argument("--umi_tag", type=str, default="UB")
    sc_g.add_argument("--strip_barcode_suffix", action="store_true", default=False)
    sc_g.add_argument("--barcode2barcode", type=str, default=None)

    algo_g.add_argument("--report_novel_unspliced", "-u", type=_bool_str, default=None)
    algo_g.add_argument(
        "--report_canonical",
        type=str,
        choices=[e.name for e in StrandnessReportingLevel],
        default=StrandnessReportingLevel.auto.name,
    )
    algo_g.add_argument(
        "--polya_requirement",
        type=str,
        choices=[e.name for e in PolyAUsageStrategies],
        default=PolyAUsageStrategies.auto.name,
    )
    algo_g.add_argument("--novel_apa", action="store_true", default=False)
    algo_g.add_argument("--no_splice_site_correction", action="store_true", default=False)
    algo_g.add_argument(
        "--transcript_quantification",
        choices=COUNTING_STRATEGIES,
        type=str,
        default=CountingStrategy_default(),
    )
    algo_g.add_argument(
        "--gene_quantification",
        choices=COUNTING_STRATEGIES,
        type=str,
        default=GeneStrategy_default(),
    )
    algo_g.add_argument(
        "--matching_strategy", choices=["exact", "precise", "default", "loose"], type=str, default=None
    )
    algo_g.add_argument(
        "--splice_correction_strategy",
        choices=["none", "default_pacbio", "default_ont", "conservative_ont", "all", "assembly"],
        type=str,
        default=None,
    )
    algo_g.add_argument(
        "--model_construction_strategy",
        choices=[
            "reliable", "default_pacbio", "sensitive_pacbio", "pacbio_all",
            "fl_pacbio", "default_ont", "sensitive_ont", "all", "assembly",
        ],
        type=str,
        default=None,
    )
    algo_g.add_argument("--delta", type=int, default=None)
    algo_g.add_argument("--use_replicas", type=_bool_str, default=True)

    parser.add_argument("--graph_clustering_distance", type=int, default=None)
    parser.add_argument("--cage", type=str, default=None)
    parser.add_argument("--cage-shift", type=int, default=50)
    parser.add_argument("--collect_polya_training", type=str, default=None)
    parser.add_argument("--collect_tss_training", type=str, default=None)
    return parser


def CountingStrategy_default() -> str:
    from isoquant_next.config.enums import CountingStrategy

    return CountingStrategy.unique_only.name


def GeneStrategy_default() -> str:
    from isoquant_next.config.enums import CountingStrategy

    return CountingStrategy.unique_splicing_consistent.name


def _validate_inputs(ns: argparse.Namespace) -> int:
    if not ns.reference or not os.path.isfile(ns.reference):
        logger.critical("Reference genome not provided or missing: %s", ns.reference)
        return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND if ns.reference else IsoQuantExitCode.INVALID_PARAMETER)
    if not ns.data_type or ns.data_type not in DATA_TYPE_ALIASES:
        logger.critical("Data type not provided or unsupported: %s", ns.data_type)
        return int(IsoQuantExitCode.INVALID_PARAMETER)
    if not any([ns.fastq, ns.bam, ns.unmapped_bam, ns.read_assignments, ns.yaml]):
        logger.critical("No input data provided")
        return int(IsoQuantExitCode.NO_INPUT_DATA)
    if ns.yaml and ns.illumina_bam:
        logger.critical("--yaml cannot be combined with --illumina_bam")
        return int(IsoQuantExitCode.INVALID_PARAMETER)
    if ns.read_group and "none" in ns.read_group and len(ns.read_group) > 1:
        logger.critical("--read_group 'none' cannot be combined with other values")
        return int(IsoQuantExitCode.INVALID_PARAMETER)
    if ns.large_output and "none" in ns.large_output and len(ns.large_output) > 1:
        logger.critical("--large_output 'none' cannot be combined with other values")
        return int(IsoQuantExitCode.INVALID_PARAMETER)

    if ns.genedb is not None and not os.path.isfile(ns.genedb):
        logger.critical("Gene database does not exist: %s", ns.genedb)
        return int(IsoQuantExitCode.GENE_DB_NOT_FOUND)
    if ns.cage is not None and not os.path.isfile(ns.cage):
        logger.critical("CAGE BED file does not exist: %s", ns.cage)
        return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
    if ns.junc_bed_file is not None and not os.path.isfile(ns.junc_bed_file):
        logger.critical("Junction BED file does not exist: %s", ns.junc_bed_file)
        return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)

    for grp_file_list in (ns.barcoded_reads, ns.barcode_whitelist):
        if grp_file_list:
            for fp in grp_file_list:
                if not os.path.isfile(fp):
                    logger.critical("Missing barcode file: %s", fp)
                    return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
    if ns.molecule and not os.path.isfile(ns.molecule):
        logger.critical("Missing molecule file: %s", ns.molecule)
        return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
    if ns.barcode2spot:
        bp, _, _ = parse_barcode2spot_spec(ns.barcode2spot)
        if not os.path.isfile(bp):
            logger.critical("Missing barcode2spot file: %s", bp)
            return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
    if ns.barcode2barcode:
        bp, _, _ = parse_barcode2spot_spec(ns.barcode2barcode)
        if not os.path.isfile(bp):
            logger.critical("Missing barcode2barcode file: %s", bp)
            return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
    if ns.read_group:
        for spec in ns.read_group:
            parts = spec.split(":")
            if len(parts) >= 2 and parts[0] == "file" and not os.path.isfile(parts[1]):
                logger.critical("Missing read_group file: %s", parts[1])
                return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)

    for lst_opt, is_bam in ((ns.fastq, False), (ns.bam, True), (ns.unmapped_bam, False), (ns.illumina_bam, True)):
        if not lst_opt:
            continue
        for fp in lst_opt:
            if not os.path.isfile(fp):
                logger.critical("Input file does not exist: %s", fp)
                return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
            if is_bam:
                with pysam.AlignmentFile(fp, "rb") as af:
                    if not af.has_index():
                        logger.critical("BAM file %s is not indexed", fp)
                        return int(IsoQuantExitCode.BAM_NOT_INDEXED)

    if ns.yaml and not os.path.isfile(ns.yaml):
        logger.critical("YAML file does not exist: %s", ns.yaml)
        return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
    if ns.read_assignments:
        for pref in ns.read_assignments:
            if not glob.glob(pref + "*"):
                logger.critical("No assignment save files with prefix %s", pref)
                return int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
    return 0


def _detect_umi_length(ns: argparse.Namespace, mode: IsoQuantMode) -> int:
    if not mode.needs_barcode_calling():
        return 0
    bc_sources = sum([bool(ns.barcode_whitelist), bool(ns.barcoded_reads), bool(ns.barcoded_bam)])
    if bc_sources > 1:
        raise IsoQuantError(
            IsoQuantExitCode.INVALID_PARAMETER,
            "--barcode_whitelist, --barcoded_reads, and --barcoded_bam are mutually exclusive",
        )
    if mode == IsoQuantMode.custom_sc and not any([ns.molecule, ns.barcoded_reads, ns.barcoded_bam]):
        raise IsoQuantError(IsoQuantExitCode.BARCODE_WHITELIST_MISSING, "custom_sc requires --molecule/--barcoded_reads/--barcoded_bam")
    if mode != IsoQuantMode.custom_sc and not any([ns.barcode_whitelist, ns.barcoded_reads, ns.barcoded_bam]):
        raise IsoQuantError(IsoQuantExitCode.BARCODE_WHITELIST_MISSING, f"Mode {mode.name} requires barcode whitelist/reads/BAM tags")

    if ns.barcoded_bam:
        bams = ns.bam or []
        if ns.yaml:
            storage = InputDataStorage(ns)
            bams = [storage.samples[0].file_list[0][0]]
        if bams and os.path.isfile(bams[0]):
            with pysam.AlignmentFile(bams[0], "rb") as af:
                for aln in af:
                    if aln.has_tag(ns.umi_tag):
                        return len(str(aln.get_tag(ns.umi_tag)))
        return 0
    from isoquant_next.barcodes.detect_barcodes import get_umi_length

    return get_umi_length(mode)


def _normalize_paths_in_ns(ns: argparse.Namespace) -> None:
    for opt in ("reference", "genedb", "index", "junc_bed_file", "cage", "genedb_output", "yaml", "molecule"):
        setattr(ns, opt, _abspath_opt(getattr(ns, opt)))
    for opt in ("bam", "fastq", "unmapped_bam", "illumina_bam", "barcode_whitelist", "barcoded_reads", "read_assignments"):
        setattr(ns, opt, _abspath_list(getattr(ns, opt)))
    ns.output = os.path.abspath(ns.output)
    if ns.read_group:
        updated = []
        for spec in ns.read_group:
            parts = spec.split(":")
            if len(parts) >= 2 and parts[0] == "file":
                parts[1] = os.path.abspath(parts[1])
            updated.append(":".join(parts))
        ns.read_group = list(dict.fromkeys(updated))


def _save_cli_options_atomic(path: Path, payload: Dict[str, Any]) -> None:
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _resolve_analyses(ns: argparse.Namespace, mode: IsoQuantMode) -> tuple[AnalysisType, ...]:
    if ns.analysis is not None:
        analyses = {ANALYSIS_ALIASES[a] for a in ns.analysis}
    elif mode.needs_barcode_calling():
        analyses = {AnalysisType.quantification}
    elif ns.genedb:
        analyses = {AnalysisType.quantification, AnalysisType.transcript_discovery}
    else:
        analyses = {AnalysisType.transcript_discovery}

    if ns.count_exons or ns.count_intron_retentions:
        analyses.add(AnalysisType.exon_quantification)
    if ns.no_model_construction:
        analyses.discard(AnalysisType.transcript_discovery)
    if not ns.genedb:
        for dep in (AnalysisType.quantification, AnalysisType.exon_quantification, AnalysisType.fusion):
            analyses.discard(dep)
        if not analyses and not ns.no_model_construction:
            analyses.add(AnalysisType.transcript_discovery)
    return tuple(sorted(analyses, key=lambda a: a.name))


def _build_params_from_ns(ns: argparse.Namespace, cmd_line: str) -> tuple[PipelineConfig, PipelineParamsView]:
    tmp_storage = InputDataStorage(ns)
    has_replicas = tmp_storage.has_replicas()
    mode = IsoQuantMode[ns.mode] if not isinstance(ns.mode, IsoQuantMode) else ns.mode
    analyses = _resolve_analyses(ns, mode)

    cfg = PipelineConfig.resolve(
        data_type=ns.data_type,
        mode=mode,
        fl_data=ns.fl_data,
        reference=ns.reference,
        genedb=ns.genedb,
        complete_genedb=ns.complete_genedb,
        matching_strategy=ns.matching_strategy,
        delta=ns.delta,
        splice_correction_strategy=ns.splice_correction_strategy,
        model_construction_strategy=ns.model_construction_strategy,
        graph_clustering_distance=ns.graph_clustering_distance,
        report_novel_unspliced=ns.report_novel_unspliced,
        polya_requirement=ns.polya_requirement,
        polya_trimmed=ns.polya_trimmed,
        report_canonical=ns.report_canonical,
        read_group=ns.read_group,
        has_replicas=has_replicas,
        barcode2spot=ns.barcode2spot,
        barcode2barcode=ns.barcode2barcode,
        bam_tags=ns.bam_tags,
        large_output=ns.large_output,
        analyses=analyses,
        threads=ns.threads,
        output_dir=ns.output,
        prefix=ns.prefix,
        resume=ns.resume,
        clean_start=ns.clean_start,
    )
    cfg = replace(
        cfg,
        quantification=QuantificationConfig(
            transcript_strategy=CountingStrategy[ns.transcript_quantification],
            gene_strategy=CountingStrategy[ns.gene_quantification],
            count_exons=AnalysisType.exon_quantification in analyses,
            count_intron_retentions=AnalysisType.exon_quantification in analyses,
        ),
    )
    umi_len = _detect_umi_length(ns, cfg.barcodes.mode)
    params = PipelineParamsView.from_config(
        cfg,
        cmd_line=cmd_line,
        fastq=ns.fastq,
        bam=ns.bam,
        unmapped_bam=ns.unmapped_bam,
        yaml=ns.yaml,
        illumina_bam=ns.illumina_bam,
        read_assignments=ns.read_assignments,
        labels=ns.labels,
        stranded=ns.stranded,
        barcode_whitelist=ns.barcode_whitelist,
        barcoded_reads=ns.barcoded_reads,
        barcoded_bam=ns.barcoded_bam,
        barcode_tag=ns.barcode_tag,
        umi_tag=ns.umi_tag,
        strip_barcode_suffix=ns.strip_barcode_suffix,
        molecule=ns.molecule,
        umi_length=umi_len,
        sqanti_output=ns.sqanti_output,
        check_canonical=ns.check_canonical,
        counts_format=ns.counts_format,
        normalization_method=ns.normalization_method,
        emit_read_ids=ns.emit_read_ids,
        old_exon_count_format=ns.old_exon_count_format,
        use_replicas=ns.use_replicas,
        cage=ns.cage,
        cage_shift=ns.cage_shift,
        process_only_chr=ns.process_only_chr,
        discard_chr=ns.discard_chr,
        keep_tmp=ns.keep_tmp,
        run_aligner_only=ns.run_aligner_only,
        gtf_check=ns.gtf_check,
        high_memory=ns.high_memory,
        gzipped=ns.gzipped,
        use_secondary=ns.use_secondary,
        min_mapq=ns.min_mapq,
        inconsistent_mapq_cutoff=ns.inconsistent_mapq_cutoff,
        simple_alignments_mapq_cutoff=ns.simple_alignments_mapq_cutoff,
        max_coverage_small_chr=ns.max_coverage_small_chr,
        max_coverage_normal_chr=ns.max_coverage_normal_chr,
        no_junc_bed=ns.no_junc_bed,
        junc_bed_file=ns.junc_bed_file,
        index=ns.index,
        aligner=ns.aligner,
        indexing_options=ns.indexing_options,
        mapping_options=ns.mapping_options,
        genedb_output=ns.genedb_output,
        collect_polya_training=ns.collect_polya_training,
        collect_tss_training=ns.collect_tss_training,
        debug=ns.debug,
        force=ns.force,
        novel_apa=ns.novel_apa,
        no_splice_site_correction=ns.no_splice_site_correction,
    )
    return cfg, params


def main(cmd_args: Optional[List[str]] = None) -> int:
    pin_single_thread_math_env(force=True)
    raw_argv = list(sys.argv[1:] if cmd_args is None else cmd_args)
    parser = build_parser()
    if not raw_argv:
        parser.print_help()
        return 0

    ns = parser.parse_args(raw_argv)
    if ns.test:
        simple_dir = Path(__file__).resolve().parents[4] / "isoquant_tests" / "simple_data"
        out_dir = os.path.abspath("isoquant_next_test")
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        return main([
            "--clean_start",
            "--data_type", "nanopore",
            "--fastq", str(simple_dir / "chr9.4M.ont.sim.fq.gz"),
            "--genedb", str(simple_dir / "chr9.4M.gtf.gz"),
            "--complete_genedb",
            "-r", str(simple_dir / "chr9.4M.fa.gz"),
            "-o", out_dir,
            "-t", "2",
        ])

    out_dir = Path(os.path.abspath(ns.output))
    manifest_path = out_dir / _MANIFEST_FILE
    cli_opts_path = out_dir / _CLI_OPTIONS_FILE

    if ns.resume:
        if not manifest_path.exists() or not cli_opts_path.exists():
            configure_logging(debug=ns.debug)
            logger.critical("Previous run manifest not found in %s; cannot --resume", out_dir)
            return int(IsoQuantExitCode.RESUME_CONFIG_NOT_FOUND)
        cfg = PipelineConfig.load_manifest(manifest_path)
        saved_ns_dict = json.loads(cli_opts_path.read_text(encoding="utf-8"))
        saved_ns_dict["resume"] = True
        saved_ns_dict["clean_start"] = False
        saved_ns_dict["threads"] = ns.threads
        saved_ns_dict["debug"] = ns.debug
        saved_ns_dict["high_memory"] = ns.high_memory
        saved_ns_dict["output"] = str(out_dir)
        ns = argparse.Namespace(**saved_ns_dict)
    else:
        _normalize_paths_in_ns(ns)

    out_dir.mkdir(parents=True, exist_ok=True)
    configure_logging(debug=ns.debug, log_file=out_dir / "isoquant.log")

    try:
        rc = _validate_inputs(ns)
        if rc != 0:
            return rc
        cmd_line = " ".join(["isoquant-next"] + raw_argv)
        cfg, params = _build_params_from_ns(ns, cmd_line=cmd_line)
        if not ns.resume:
            cfg.save_manifest(manifest_path)
            _save_cli_options_atomic(cli_opts_path, vars(ns))
        run_pipeline(params, config_fingerprint=cfg.command_fingerprint())
        return 0
    except IsoQuantError as err:
        logger.critical("%s", err)
        return int(err.code)


if __name__ == "__main__":
    raise SystemExit(main())
