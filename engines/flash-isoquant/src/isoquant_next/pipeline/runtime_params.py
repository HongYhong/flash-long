"""Picklable runtime parameter bridge, `.iqidx` string pools builder, and BED12 extractor."""

from __future__ import annotations

from dataclasses import asdict
import glob
import os
from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Optional, Tuple

from isoquant_next import __version__
from isoquant_next.alignment.read_groups import (
    get_grouping_pool_types,
    parse_barcode2spot_spec,
)
from isoquant_next.annotation.iqidx import MmapAnnotationIndex
from isoquant_next.config.enums import (
    ASSEMBLY,
    AnalysisType,
    NANOPORE_DATA,
    PACBIO_CCS_DATA,
    PolyATrimmed,
    PolyAUsageStrategies,
    StrandnessReportingLevel,
    set_polya_requirement_strategy,
)
from isoquant_next.config.settings import PipelineConfig
from isoquant_next.core.error_codes import IsoQuantError, IsoQuantExitCode
from isoquant_next.core.string_pools import StringPoolManager
from isoquant_next.pipeline.input_data import InputDataStorage, InputDataType


class PipelineParamsView:
    """Mutable, shallow-copyable parameter view hydrated from PipelineConfig and CLI options."""

    def __init__(self, **kwargs: Any) -> None:
        self.__dict__.update(kwargs)

    @classmethod
    def from_config(
        cls,
        cfg: PipelineConfig,
        *,
        cmd_line: str = "",
        fastq: Optional[List[str]] = None,
        bam: Optional[List[str]] = None,
        unmapped_bam: Optional[List[str]] = None,
        yaml: Optional[str] = None,
        illumina_bam: Optional[List[str]] = None,
        read_assignments: Optional[List[str]] = None,
        labels: Optional[List[str]] = None,
        stranded: str = "none",
        barcode_whitelist: Optional[List[str]] = None,
        barcoded_reads: Optional[List[str]] = None,
        barcoded_bam: bool = False,
        barcode_tag: str = "CB",
        umi_tag: str = "UB",
        strip_barcode_suffix: bool = False,
        molecule: Optional[str] = None,
        umi_length: int = 0,
        sqanti_output: bool = False,
        check_canonical: bool = False,
        counts_format: Optional[List[str]] = None,
        normalization_method: str = "simple",
        emit_read_ids: bool = False,
        old_exon_count_format: bool = False,
        use_replicas: bool = True,
        cage: Optional[str] = None,
        cage_shift: int = 50,
        process_only_chr: Optional[List[str]] = None,
        discard_chr: Optional[List[str]] = None,
        keep_tmp: bool = False,
        run_aligner_only: bool = False,
        gtf_check: bool = True,
        high_memory: bool = False,
        gzipped: bool = True,
        use_secondary: bool = False,
        min_mapq: Optional[int] = None,
        inconsistent_mapq_cutoff: int = 5,
        simple_alignments_mapq_cutoff: int = 1,
        max_coverage_small_chr: int = 1_000_000,
        max_coverage_normal_chr: int = -1,
        no_junc_bed: bool = False,
        junc_bed_file: Optional[str] = None,
        index: Optional[str] = None,
        aligner: Optional[str] = None,
        indexing_options: Optional[str] = None,
        mapping_options: Optional[str] = None,
        genedb_output: Optional[str] = None,
        collect_polya_training: Optional[str] = None,
        collect_tss_training: Optional[str] = None,
        debug: bool = False,
        force: bool = False,
        novel_apa: bool = False,
        no_splice_site_correction: bool = False,
        umi_dedup_min_count: int = 1,
    ) -> "PipelineParamsView":
        m = cfg.matching
        sc = cfg.splice_correction
        mc = cfg.model_construction
        q = cfg.quantification
        bc = cfg.barcodes
        rt = cfg.runtime

        analyses = set(cfg.analyses)
        run_quantification = AnalysisType.quantification in analyses
        count_exons = AnalysisType.exon_quantification in analyses
        count_intron_retentions = count_exons
        fusion = AnalysisType.fusion in analyses
        no_model_construction = AnalysisType.transcript_discovery not in analyses
        predict_terminal_sites = bool(cfg.genedb) and (run_quantification or not no_model_construction)

        is_ont = cfg.data_type == NANOPORE_DATA
        needs_reference = True

        flat_cfg = {**asdict(m), **asdict(sc), **asdict(mc)}
        flat_cfg.update(
            dict(
                matching_strategy_name=m.strategy_name,
                splice_correction_strategy_name=sc.strategy_name,
                model_construction_strategy_name=mc.strategy_name,
                data_type=cfg.data_type,
                fl_data=cfg.fl_data,
                fl_only=mc.fl_only,
                reference=cfg.reference,
                genedb=cfg.genedb,
                original_annotation=cfg.genedb,
                complete_genedb=cfg.complete_genedb,
                analyses=analyses,
                run_quantification=run_quantification,
                count_exons=count_exons,
                count_intron_retentions=count_intron_retentions,
                fusion=fusion,
                no_model_construction=no_model_construction,
                predict_terminal_sites=predict_terminal_sites,
                transcript_quantification=q.transcript_strategy.name,
                gene_quantification=q.gene_strategy.name,
                counts_format=list(counts_format or ["default"]),
                normalization_method=normalization_method,
                emit_read_ids=emit_read_ids,
                old_exon_count_format=old_exon_count_format,
                requires_polya_for_construction=set_polya_requirement_strategy(
                    False, mc.polya_requirement_strategy
                ),
                require_monointronic_polya=mc.require_monointronic_polya,
                require_monoexonic_polya=mc.require_monoexonic_polya,
                report_canonical_strategy=mc.report_canonical,
                polya_requirement_strategy=mc.polya_requirement_strategy,
                cage=cage,
                cage_shift=cage_shift,
                novel_apa=novel_apa,
                no_splice_site_correction=no_splice_site_correction,
                correct_splice_sites=not no_splice_site_correction,
                mode=bc.mode,
                barcoded_reads=barcoded_reads,
                barcode_whitelist=barcode_whitelist,
                barcoded_bam=barcoded_bam,
                barcode_tag=barcode_tag,
                umi_tag=umi_tag,
                strip_barcode_suffix=strip_barcode_suffix,
                molecule=molecule,
                umi_length=umi_length,
                barcode2spot=bc.barcode2spot,
                barcode2barcode=bc.barcode2barcode,
                read_group=list(bc.read_group),
                large_output=list(bc.large_output) if bc.large_output else ["none"],
                bam_tags=list(bc.bam_tags),
                fastq=fastq,
                bam=bam,
                unmapped_bam=unmapped_bam,
                yaml=yaml,
                illumina_bam=illumina_bam,
                read_assignments=read_assignments,
                labels=labels,
                stranded=stranded,
                output=rt.output_dir,
                prefix=rt.prefix,
                threads=rt.threads,
                clean_start=rt.clean_start,
                resume=rt.resume,
                force=force,
                debug=debug,
                genedb_output=genedb_output or rt.output_dir,
                keep_tmp=keep_tmp,
                run_aligner_only=run_aligner_only,
                gtf_check=gtf_check,
                high_memory=high_memory,
                gzipped=gzipped,
                needs_reference=needs_reference,
                use_secondary=use_secondary,
                min_mapq=min_mapq,
                inconsistent_mapq_cutoff=inconsistent_mapq_cutoff,
                simple_alignments_mapq_cutoff=simple_alignments_mapq_cutoff,
                max_coverage_small_chr=max_coverage_small_chr,
                max_coverage_normal_chr=max_coverage_normal_chr,
                no_junc_bed=no_junc_bed or (cfg.genedb is None),
                junc_bed_file=junc_bed_file,
                index=index,
                aligner=aligner,
                indexing_options=indexing_options,
                mapping_options=mapping_options,
                use_replicas=use_replicas,
                sqanti_output=sqanti_output and bool(cfg.genedb) and not no_model_construction,
                check_canonical=check_canonical,
                collect_polya_training=collect_polya_training,
                collect_tss_training=collect_tss_training,
                umi_dedup_min_count=umi_dedup_min_count,
                fai_file_name=f"{cfg.reference}.fai" if cfg.reference else None,
                alignment_config_path=os.path.join(rt.output_dir, "alignment_config.json"),
                log_file=os.path.join(rt.output_dir, "isoquant.log"),
                _cmd_line=cmd_line,
                _version=__version__,
            )
        )
        view = cls(**flat_cfg)
        # Ensure minor_exon_extension / spliced-read defaults match isoquant.py
        if not cfg.genedb and not is_ont and not hasattr(view, "correct_minor_errors"):
            view.correct_minor_errors = False
        view.input_data = InputDataStorage(view)
        return view


def setup_string_pools(
    args,
    sample,
    chr_ids: Iterable[str],
    chr_id: Optional[str] = None,
    genedb_index=None,
    *,
    load_barcode_pool: bool = False,
    load_tsv_pools: bool = False,
    read_group_file_prefix: Optional[str] = None,
) -> StringPoolManager:
    """Construct StringPoolManager directly from mmap `.iqidx` (zero SQLite overhead)."""
    sp = StringPoolManager()
    sp.build_chromosome_pool(chr_ids)

    close_db = False
    db = genedb_index
    if db is None and getattr(args, "genedb", None):
        if str(args.genedb).endswith(".iqidx"):
            db = MmapAnnotationIndex(args.genedb)
            close_db = True

    try:
        if isinstance(db, MmapAnnotationIndex):
            sp.build_from_iqidx(db)
        elif db is not None:
            sp.build_from_gffutils(db)
    finally:
        if close_db and db is not None:
            db.close()

    pool_types = get_grouping_pool_types(args)
    for spec_idx, pool_type in pool_types.items():
        sp.set_group_spec_pool_type(spec_idx, pool_type)

    if any(pt == "file_name" for pt in pool_types.values()):
        sp.build_file_name_pool(sample)

    if any(pt == "barcode_spot" or pt.startswith("barcode_spot:") for pt in pool_types.values()):
        if getattr(args, "barcode2spot", None):
            _fn, _bcol, spot_cols = parse_barcode2spot_spec(args.barcode2spot)
            for col_idx, file_col in enumerate(spot_cols):
                sp.build_barcode_spot_pool(args.barcode2spot, column_index=col_idx, file_column=file_col)

    if chr_id is not None:
        if load_barcode_pool and getattr(sample, "barcodes_split_reads", None):
            barcode_file = sample.get_barcodes_split_file(chr_id)
            sp.load_barcode_pool(barcode_file)
            if len(sp.barcode_pool) == 0 and getattr(args, "barcoded_bam", False):
                sp.load_barcode_umi_pools(sample.get_barcode_pools_file(chr_id))

        rg_file = read_group_file_prefix if read_group_file_prefix else getattr(sample, "read_group_file", None)
        if load_tsv_pools and rg_file:
            tsv_pool_info: Dict[int, List[Tuple[int, str, str]]] = {}
            for pt in pool_types.values():
                if pt.startswith("tsv:"):
                    parts = pt.split(":")
                    if len(parts) >= 4:
                        tsv_spec = int(parts[1])
                        col_idx = int(parts[2])
                        delim = parts[3]
                        pool_key = f"{tsv_spec}:{col_idx}"
                        tsv_pool_info.setdefault(tsv_spec, []).append((col_idx, delim, pool_key))

            pattern = f"{rg_file}_spec*_{chr_id}"
            for spec_file in sorted(glob.glob(pattern)):
                m = re.search(r"_spec(\d+)_", spec_file)
                if m:
                    spec_idx = int(m.group(1))
                    for col_idx, delim, pool_key in tsv_pool_info.get(spec_idx, []):
                        sp.load_read_group_tsv_pool(spec_file, pool_key, col_index=col_idx, delimiter=delim)

    return sp


def iqidx_to_bed12(iqidx_path: str | Path, out_bed: str | Path) -> Path:
    """Convert zero-copy `.iqidx` annotation into BED12 junction file for minimap2 `--junc-bed`."""
    from isoquant_next.annotation.gtf_validator import parse_gff3_attributes, parse_gtf_attributes

    _colors = {
        "protein_coding": "0,128,255",
        "mRNA": "0,128,255",
        "lincRNA": "0,192,0",
        "snRNA": "0,192,0",
        "miRNA": "0,192,0",
        "misc_RNA": "0,192,0",
    }
    _strand = {1: "+", -1: "-", 0: "."}

    out_p = Path(out_bed)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    tmp_p = out_p.with_suffix(out_p.suffix + f".tmp.{os.getpid()}")

    with MmapAnnotationIndex(iqidx_path) as idx, open(tmp_p, "w", encoding="utf-8") as fh:
        g_tx = idx.gene_tx_ptr
        g_id_s = idx.gene_id_str
        g_attr_s = idx.gene_attr_str
        t_id_s = idx.tx_id_str
        t_str = idx.tx_strand
        t_attr_s = idx.tx_attr_str
        t_ex = idx.tx_exon_ptr
        t_ot = idx.tx_other_ptr
        e_start = idx.exon_start
        e_end = idx.exon_end
        o_start = idx.other_start
        o_end = idx.other_end
        o_type = idx.other_type

        def _parse_attr(raw_s: str) -> Dict[str, str]:
            if not raw_s:
                return {}
            return parse_gff3_attributes(raw_s) if ("=" in raw_s and '"' not in raw_s) else parse_gtf_attributes(raw_s)

        for chrom in idx.seqids():
            g0, g1 = idx.chrom_gene_slice(chrom)
            for gi in range(g0, g1):
                gid = idx.get_str(int(g_id_s[gi]))
                gattrs = _parse_attr(idx.get_str(int(g_attr_s[gi])))
                default_gname = gattrs.get("gene_name") or gid

                for ti in range(int(g_tx[gi]), int(g_tx[gi + 1])):
                    ex0, ex1 = int(t_ex[ti]), int(t_ex[ti + 1])
                    if ex0 >= ex1:
                        continue
                    tid = idx.get_str(int(t_id_s[ti]))
                    tattrs = _parse_attr(idx.get_str(int(t_attr_s[ti])))
                    ttype = tattrs.get("transcript_type", "")
                    gname = tattrs.get("gene_name") or default_gname
                    tname = f"{tid}|{ttype}|{gname}"
                    strand = _strand.get(int(t_str[ti]), ".")
                    exons = [(int(e_start[ei]) - 1, int(e_end[ei])) for ei in range(ex0, ex1)]
                    tx_start = exons[0][0]
                    tx_end = exons[-1][1]

                    cds_coords = [
                        (int(o_start[oi]) - 1, int(o_end[oi]))
                        for oi in range(int(t_ot[ti]), int(t_ot[ti + 1]))
                        if int(o_type[oi]) == 0
                    ]
                    if cds_coords:
                        thick_start, thick_end = cds_coords[0][0], cds_coords[-1][1]
                    else:
                        thick_start, thick_end = tx_start, tx_end

                    blens = ",".join(str(e - s) for s, e in exons) + ","
                    bstarts = ",".join(str(s - tx_start) for s, _e in exons) + ","
                    color = _colors.get(ttype, "196,196,196")
                    fh.write(
                        f"{chrom}\t{tx_start}\t{tx_end}\t{tname}\t1000\t{strand}\t"
                        f"{thick_start}\t{thick_end}\t{color}\t{len(exons)}\t{blens}\t{bstarts}\n"
                    )
    os.replace(tmp_p, out_p)
    return out_p
