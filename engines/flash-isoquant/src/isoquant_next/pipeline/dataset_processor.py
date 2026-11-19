"""O(Locus) streaming & two-stage DatasetProcessor with zero cross-process ReadAssignment pickling."""

from __future__ import annotations

from collections import defaultdict
import gc
import glob
import gzip
import logging
import os
import shutil
from typing import Dict, List, Optional, Tuple

import pysam

from isoquant_next.alignment.read_groups import (
    create_read_grouper,
    get_grouping_pool_types,
    get_grouping_strategy_names,
    load_barcode2barcode_mapping,
    parse_barcode2spot_spec,
    prepare_read_groups,
)
from isoquant_next.annotation.gene_info import TranscriptModelType
from isoquant_next.annotation.iqidx import MmapAnnotationIndex
from isoquant_next.annotation.reference import ReferenceGenome
from isoquant_next.assignment.io import IOSupport, ReadInfoPrinter, SqantiTSVPrinter, TmpFileAssignmentPrinter, VoidPrinter
from isoquant_next.assignment.loader import create_assignment_loader
from isoquant_next.assignment.multimap import MultimapResolver, MultimapResolvingStrategy
from isoquant_next.assignment.types import BasicReadAssignment, CompactMultimapTuple, ReadAssignmentType
from isoquant_next.barcodes.table_splitter import split_read_table_by_chromosomes
from isoquant_next.barcodes.umi_filtering import UMIFilter, create_transcript_info_dict
from isoquant_next.config.enums import IsoQuantMode, PolyAUsageStrategies, set_polya_requirement_strategy
from isoquant_next.core.concurrency import WorkerPoolExecutor
from isoquant_next.core.error_codes import IsoQuantError, IsoQuantExitCode
from isoquant_next.core.file_naming import (
    allinfo_file_name,
    allinfo_stats_file_name,
    bamstat_file_name,
    convert_chr_id_to_file_name_str,
    dynamic_pools_file_name,
    info_file_name,
    read_group_file_from_saves,
    read_group_lock_filename,
    read_groups_file_name,
    read_stat_file_name,
    reads_collected_lock_file_name,
    reads_processed_lock_file_name,
    saves_file_name,
    split_barcodes_lock_filename,
    transcript_stat_file_name,
    umi_barcode2barcode_global_lock,
    umi_barcode2barcode_prefix,
    umi_filtered_global_lock_file_name,
    umi_filtered_lock_file_name,
    umi_output_prefix,
)
from isoquant_next.core.file_utils import merge_counts, merge_files
from isoquant_next.core.id_policy import ExcludingIdDistributor, FeatureIdStorage, SimpleIDDistributor
from isoquant_next.core.serialization import read_int, read_list, read_string, write_int, write_list, write_string
from isoquant_next.core.stats import EnumStats
from isoquant_next.core.string_pools import StringPoolManager
from isoquant_next.core.utils import large_output_enabled, proper_plural_form
from isoquant_next.discovery.builder import GraphBasedModelConstructor
from isoquant_next.discovery.printer import GFFPrinter, VoidTranscriptPrinter, create_extended_storage
from isoquant_next.discovery.state import ChromosomeDiscoveryState
from isoquant_next.pipeline.aggregator import ReadAssignmentAggregator
from isoquant_next.pipeline.alignment_collector import AlignmentCollector, AlignmentType
from isoquant_next.pipeline.input_data import SampleData
from isoquant_next.pipeline.multimap_resolver import (
    prepare_compact_multimappers_from_saves,
    resolve_compact_multimapper_records,
)
from isoquant_next.pipeline.runtime_params import setup_string_pools

logger = logging.getLogger("IsoQuant")

_UMI_ED_BY_MODE = {
    IsoQuantMode.bulk: [],
    IsoQuantMode.tenX_v3: [3],
    IsoQuantMode.tenX_v2: [3],
    IsoQuantMode.tenX_v3_split: [3],
    IsoQuantMode.tenX_v2_split: [3],
    IsoQuantMode.visium_5prime: [3],
    IsoQuantMode.curio: [3],
    IsoQuantMode.visium_hd: [4],
    IsoQuantMode.stereoseq: [4],
    IsoQuantMode.stereoseq_nosplit: [4],
    IsoQuantMode.custom_sc: [4],
}


def _save_dynamic_pools(sp: StringPoolManager, path: str) -> None:
    if sp.has_dynamic_pools():
        with open(path, "wb") as fh:
            sp.serialize_dynamic_pools(fh)


def _load_dynamic_pools(sp: StringPoolManager, path: str) -> None:
    if os.path.exists(path):
        with open(path, "rb") as fh:
            sp.deserialize_dynamic_pools(fh)


def _save_read_groups(path: str, read_groups) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        if isinstance(read_groups, list):
            fh.write(f"{len(read_groups)}\n")
            for gs in read_groups:
                fh.write(";".join(sorted(gs)) + "\n")
        else:
            fh.write("1\n")
            fh.write(";".join(sorted(read_groups)) + "\n")


def _load_read_groups(path: str, read_grouper) -> None:
    with open(path, "r", encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh]
    n = int(lines[0])
    if isinstance(read_grouper.read_groups, list):
        read_grouper.read_groups = [set(lines[i + 1].split(";")) if lines[i + 1] else set() for i in range(n)]
    else:
        gstr = lines[1] if len(lines) > 1 else ""
        read_grouper.read_groups = set(gstr.split(";")) if gstr else set()


def _load_barcode_dict(sample: SampleData, chr_id: str) -> Dict[str, Tuple[str, str]]:
    out: Dict[str, Tuple[str, str]] = {}
    if not sample.barcodes_split_reads:
        return out
    bf = sample.get_barcodes_split_file(chr_id)
    if not os.path.exists(bf):
        return out
    with open(bf, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line[0] == "#" or line.startswith("read_id\t"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) >= 3:
                out[parts[0]] = (parts[1], parts[2])
    return out


def _single_pass_chr_worker(task: Tuple[SampleData, str, List[str], object, Optional[List[Set[str]]]]):
    """Stream one chromosome locus-by-locus in O(Locus) RAM without `.save` dump/load I/O."""
    sample, chr_id, chr_ids, args, _all_read_groups_unused = task
    save_file = sample.get_save_file(chr_id)
    bamstat_file = sample.get_bamstat_file(chr_id)
    lock_file = reads_processed_lock_file_name(sample.out_raw_file, chr_id)
    chr_read_stat_file = read_stat_file_name(save_file)
    chr_transcript_stat_file = transcript_stat_file_name(save_file)
    construct_models = not args.no_model_construction

    if args.resume and os.path.exists(lock_file):
        al_stat = EnumStats()
        if os.path.exists(bamstat_file):
            al_stat.load(bamstat_file, enum_cls=AlignmentType)
        read_stat = EnumStats(chr_read_stat_file) if os.path.exists(chr_read_stat_file) else EnumStats()
        tx_stat = EnumStats(chr_transcript_stat_file) if (construct_models and os.path.exists(chr_transcript_stat_file)) else EnumStats()
        return chr_id, [], al_stat.stats_dict, read_stat.stats_dict, tx_stat.stats_dict

    if os.path.exists(lock_file):
        os.remove(lock_file)

    bam_files = [x[0] for x in sample.file_list]
    genedb = MmapAnnotationIndex(args.genedb) if args.genedb else None
    ref_dir = os.path.dirname(args.fai_file_name) if args.fai_file_name else None
    ref = ReferenceGenome(args.reference, index_dir=ref_dir)

    try:
        chr_record = ref.contig_view(chr_id)
        read_grouper = create_read_grouper(args, sample, chr_id)
        string_pools = setup_string_pools(
            args, sample, chr_ids, chr_id, genedb_index=genedb, load_barcode_pool=True, load_tsv_pools=True
        )
        barcode_dict = _load_barcode_dict(sample, chr_id)
        collector = AlignmentCollector(
            chr_id,
            bam_files,
            args,
            sample.illumina_bam,
            genedb,
            chr_record,
            read_grouper,
            barcode_dict,
            args.max_coverage_small_chr,
            args.usual_gene_max_coverage if hasattr(args, "usual_gene_max_coverage") else args.max_coverage_normal_chr,
            string_pools,
        )

        grouping_strategy_names = get_grouping_strategy_names(args)
        chr_state = ChromosomeDiscoveryState()
        io_support = IOSupport(args)
        transcript_stat_counter = EnumStats()
        chr_prefix = sample.get_chr_prefix(chr_id)
        tx_id_dist = ExcludingIdDistributor(genedb, chr_id)
        exon_id_storage = FeatureIdStorage(SimpleIDDistributor(), genedb, chr_id, "exon")
        tmp_gff = (
            GFFPrinter(sample.out_dir, chr_prefix, exon_id_storage, check_canonical=args.check_canonical)
            if construct_models
            else VoidTranscriptPrinter()
        )
        tmp_ext_gff = (
            GFFPrinter(
                sample.out_dir, chr_prefix, exon_id_storage,
                gtf_suffix=".extended_annotation.gtf", check_canonical=args.check_canonical,
            )
            if (construct_models and genedb)
            else VoidTranscriptPrinter()
        )

        with ReadAssignmentAggregator(
            args, sample, string_pools, genedb, chr_id, grouping_strategy_names=grouping_strategy_names
        ) as aggregator:
            model_reads_printer = (
                ReadInfoPrinter(
                    sample.get_transcript_model_reads_file(chr_id), args, io_support,
                    additional_header=aggregator.common_header,
                )
                if (construct_models and large_output_enabled(args, "read2transcripts"))
                else VoidPrinter()
            )
            sqanti_t2t_printer = (
                SqantiTSVPrinter(sample.get_t2t_tsv_file(chr_id), args, io_support)
                if args.sqanti_output
                else VoidTranscriptPrinter()
            )
            novel_models = []
            try:
                for gene_info, assignment_storage in collector.process():
                    for ra in assignment_storage:
                        if ra is None:
                            continue
                        aggregator.read_stat_counter.add(ra.assignment_type)
                        aggregator.global_printer.add_read_info(ra)
                        aggregator.global_counter.add_read_info(ra)
                    aggregator.global_counter.flush()

                    if construct_models:
                        polya_c = aggregator.polya_counter
                        tss_c = aggregator.tss_counter
                        mc = GraphBasedModelConstructor(
                            gene_info,
                            chr_record,
                            args,
                            aggregator.transcript_model_global_counter,
                            aggregator.gene_model_global_counter,
                            tx_id_dist,
                            grouping_strategy_names=grouping_strategy_names,
                            use_technical_replicas=sample.use_technical_replicas,
                            string_pools=string_pools,
                            polya_predictions=polya_c.last_gene_predictions if polya_c else None,
                            tss_predictions=tss_c.last_gene_predictions if tss_c else None,
                            chr_state=chr_state,
                        )
                        mc.process(assignment_storage)
                        if args.check_canonical:
                            io_support.add_canonical_info(mc.transcript_model_storage, gene_info)
                        tmp_gff.dump(mc.gene_info, mc.transcript_model_storage)
                        if large_output_enabled(args, "read2transcripts"):
                            for ma in mc.model_read_assignments:
                                model_reads_printer.add_read_info(ma)
                        for m in mc.transcript_model_storage:
                            if m.transcript_type != TranscriptModelType.known:
                                novel_models.append(m)
                            transcript_stat_counter.add(m.transcript_type)
                        for a in mc.transcript2transcript:
                            sqanti_t2t_printer.add_read_info(a)
                    # Inline locus release: drop assignment objects immediately after locus closure
                    assignment_storage.clear()

                aggregator.global_counter.dump()
                aggregator.read_stat_counter.dump(chr_read_stat_file)
                collector.alignment_stat_counter.dump(bamstat_file)
                if construct_models:
                    if genedb:
                        all_models, ext_gi = create_extended_storage(genedb, chr_id, chr_record, novel_models)
                        if args.check_canonical:
                            io_support.add_canonical_info(all_models, ext_gi)
                        tmp_ext_gff.dump(ext_gi, all_models)
                    aggregator.transcript_model_global_counter.dump()
                    aggregator.gene_model_global_counter.dump()
                    transcript_stat_counter.dump(chr_transcript_stat_file)
                    model_reads_printer.flush()
            finally:
                tmp_gff.close()
                tmp_ext_gff.close()
                model_reads_printer.close()
                sqanti_t2t_printer.close()

        open(lock_file, "w", encoding="utf-8").close()
        return (
            chr_id,
            read_grouper.read_groups,
            dict(collector.alignment_stat_counter.stats_dict),
            dict(aggregator.read_stat_counter.stats_dict),
            dict(transcript_stat_counter.stats_dict),
        )
    finally:
        ref.close()
        if genedb is not None:
            genedb.close()


def _collect_reads_worker(task: Tuple[SampleData, str, List[str], object]):
    """Two-stage worker: collect assignments into `.save_<chr>` and return only CompactMultimapTuple bytes."""
    sample, chr_id, chr_ids, args = task
    save_file = sample.get_save_file(chr_id)
    group_file = read_groups_file_name(save_file)
    bamstat_file = sample.get_bamstat_file(chr_id)
    lock_file = reads_collected_lock_file_name(sample.out_raw_file, chr_id)

    read_grouper = create_read_grouper(args, sample, chr_id)
    if args.resume and os.path.exists(lock_file):
        _load_read_groups(group_file, read_grouper)
        al_stat = EnumStats()
        al_stat.load(bamstat_file, enum_cls=AlignmentType)
        return chr_id, read_grouper.read_groups, dict(al_stat.stats_dict), True, [], [], 0, 0

    if os.path.exists(lock_file):
        os.remove(lock_file)

    bam_files = [x[0] for x in sample.file_list]
    genedb = MmapAnnotationIndex(args.genedb) if args.genedb else None
    ref_dir = os.path.dirname(args.fai_file_name) if args.fai_file_name else None
    ref = ReferenceGenome(args.reference, index_dir=ref_dir)
    compact_payloads: List[bytes] = []
    read_counts: Dict[str, int] = defaultdict(int)
    local_total = 0
    local_polya = 0

    try:
        chr_record = ref.contig_view(chr_id)
        string_pools = setup_string_pools(
            args, sample, chr_ids, chr_id, genedb_index=genedb, load_barcode_pool=True, load_tsv_pools=True
        )
        barcode_dict = _load_barcode_dict(sample, chr_id)
        collector = AlignmentCollector(
            chr_id, bam_files, args, sample.illumina_bam, genedb, chr_record, read_grouper,
            barcode_dict, args.max_coverage_small_chr, args.max_coverage_normal_chr, string_pools,
        )
        need_cross_chr_mm = bool(args.use_secondary)
        high_mem = bool(getattr(args, "high_memory", False))

        with TmpFileAssignmentPrinter(save_file, args) as tmp_printer:
            for gene_info, assignment_storage in collector.process():
                tmp_printer.add_gene_info(gene_info)
                for ra in assignment_storage:
                    tmp_printer.add_read_info(ra)
                    if not need_cross_chr_mm:
                        read_counts[ra.read_id] += 1
                    elif high_mem:
                        compact_payloads.append(BasicReadAssignment(ra).to_compact_tuple().to_bytes())
                    else:
                        read_counts[ra.read_id] += 1

        _save_read_groups(group_file, read_grouper.read_groups)
        collector.alignment_stat_counter.dump(bamstat_file)
        _save_dynamic_pools(string_pools, sample.get_dynamic_pools_file(chr_id))
        if getattr(args, "barcoded_bam", False) and len(string_pools.barcode_pool) > 0:
            string_pools.save_barcode_umi_pools(sample.get_barcode_pools_file(chr_id))

        # When secondary alignments are disabled, multimappers exist at most within this chromosome!
        if not need_cross_chr_mm:
            mm_dict, u_tot, u_pa = prepare_compact_multimappers_from_saves([chr_id], sample, read_counts, string_pools)
            mm_tot, mm_pa = resolve_compact_multimapper_records(
                [chr_id], sample, mm_dict, args.multimap_strategy, string_pools
            )
            local_total = u_tot + mm_tot
            local_polya = u_pa + mm_pa
            read_counts.clear()

        open(lock_file, "w", encoding="utf-8").close()
        return (
            chr_id,
            read_grouper.read_groups,
            dict(collector.alignment_stat_counter.stats_dict),
            False,
            compact_payloads,
            list(read_counts.items()) if need_cross_chr_mm else [],
            local_total,
            local_polya,
        )
    finally:
        ref.close()
        if genedb is not None:
            genedb.close()


def _construct_models_worker(task: Tuple[SampleData, str, List[str], str, object, List[Set[str]]]):
    sample, chr_id, chr_ids, saves_prefix, args, _all_read_groups = task
    use_filtered = args.mode.needs_pcr_deduplication()
    rg_prefix = read_group_file_from_saves(saves_prefix)
    pool_types = get_grouping_pool_types(args)
    needs_bc_pool = (
        use_filtered
        or getattr(args, "barcoded_bam", False)
        or any(
            pt == "barcode" or pt.startswith("barcode_spot") or pt.startswith("barcode_barcode")
            for pt in pool_types.values()
        )
    )
    genedb = MmapAnnotationIndex(args.genedb) if args.genedb else None
    try:
        string_pools = setup_string_pools(
            args, sample, chr_ids, chr_id, genedb_index=genedb,
            load_barcode_pool=needs_bc_pool, load_tsv_pools=True, read_group_file_prefix=rg_prefix,
        )
        _load_dynamic_pools(string_pools, dynamic_pools_file_name(saves_prefix, chr_id))

        chr_prefix = sample.get_chr_prefix(chr_id)
        chr_dump_file = saves_file_name(saves_prefix, chr_id)
        lock_file = reads_processed_lock_file_name(saves_prefix, chr_id)
        chr_read_stat_file = read_stat_file_name(chr_dump_file)
        chr_tx_stat_file = transcript_stat_file_name(chr_dump_file)
        construct_models = not args.no_model_construction

        has_chr_outputs = bool(glob.glob(os.path.join(sample.out_dir, f"{chr_prefix}.*")))
        if args.resume and os.path.exists(lock_file) and has_chr_outputs:
            rs = EnumStats(chr_read_stat_file)
            ts = EnumStats(chr_tx_stat_file) if construct_models else EnumStats()
            return dict(rs.stats_dict), dict(ts.stats_dict)
        if os.path.exists(lock_file):
            os.remove(lock_file)

        grouping_strategy_names = get_grouping_strategy_names(args)
        chr_state = ChromosomeDiscoveryState()
        tx_stat_counter = EnumStats()
        io_support = IOSupport(args)
        tx_id_dist = ExcludingIdDistributor(genedb, chr_id)
        exon_id_storage = FeatureIdStorage(SimpleIDDistributor(), genedb, chr_id, "exon")

        tmp_gff = (
            GFFPrinter(sample.out_dir, chr_prefix, exon_id_storage, check_canonical=args.check_canonical)
            if construct_models
            else VoidTranscriptPrinter()
        )
        tmp_ext_gff = (
            GFFPrinter(
                sample.out_dir, chr_prefix, exon_id_storage,
                gtf_suffix=".extended_annotation.gtf", check_canonical=args.check_canonical,
            )
            if (construct_models and genedb)
            else VoidTranscriptPrinter()
        )

        with (
            create_assignment_loader(chr_id, saves_prefix, genedb, args.reference, args.fai_file_name, string_pools, use_filtered) as loader,
            ReadAssignmentAggregator(
                args, sample, string_pools, genedb, chr_id, grouping_strategy_names=grouping_strategy_names
            ) as aggregator,
        ):
            model_reads_printer = (
                ReadInfoPrinter(
                    sample.get_transcript_model_reads_file(chr_id), args, io_support,
                    additional_header=aggregator.common_header,
                )
                if (construct_models and large_output_enabled(args, "read2transcripts"))
                else VoidPrinter()
            )
            sqanti_t2t_printer = (
                SqantiTSVPrinter(sample.get_t2t_tsv_file(chr_id), args, io_support)
                if args.sqanti_output
                else VoidTranscriptPrinter()
            )
            novel_models = []
            try:
                while loader.has_next():
                    gene_info, assignment_storage = loader.get_next()
                    for ra in assignment_storage:
                        if ra is None:
                            continue
                        aggregator.read_stat_counter.add(ra.assignment_type)
                        aggregator.global_printer.add_read_info(ra)
                        aggregator.global_counter.add_read_info(ra)
                    aggregator.global_counter.flush()

                    if construct_models:
                        polya_c = aggregator.polya_counter
                        tss_c = aggregator.tss_counter
                        mc = GraphBasedModelConstructor(
                            gene_info,
                            loader.chr_record,
                            args,
                            aggregator.transcript_model_global_counter,
                            aggregator.gene_model_global_counter,
                            tx_id_dist,
                            grouping_strategy_names=grouping_strategy_names,
                            use_technical_replicas=sample.use_technical_replicas,
                            string_pools=string_pools,
                            polya_predictions=polya_c.last_gene_predictions if polya_c else None,
                            tss_predictions=tss_c.last_gene_predictions if tss_c else None,
                            chr_state=chr_state,
                        )
                        mc.process(assignment_storage)
                        if args.check_canonical:
                            io_support.add_canonical_info(mc.transcript_model_storage, gene_info)
                        tmp_gff.dump(mc.gene_info, mc.transcript_model_storage)
                        if large_output_enabled(args, "read2transcripts"):
                            for ma in mc.model_read_assignments:
                                model_reads_printer.add_read_info(ma)
                        for m in mc.transcript_model_storage:
                            if m.transcript_type != TranscriptModelType.known:
                                novel_models.append(m)
                            tx_stat_counter.add(m.transcript_type)
                        for a in mc.transcript2transcript:
                            sqanti_t2t_printer.add_read_info(a)
                    assignment_storage.clear()

                aggregator.global_counter.dump()
                aggregator.read_stat_counter.dump(chr_read_stat_file)
                if construct_models:
                    if genedb:
                        all_models, ext_gi = create_extended_storage(genedb, chr_id, loader.chr_record, novel_models)
                        if args.check_canonical:
                            io_support.add_canonical_info(all_models, ext_gi)
                        tmp_ext_gff.dump(ext_gi, all_models)
                    aggregator.transcript_model_global_counter.dump()
                    aggregator.gene_model_global_counter.dump()
                    tx_stat_counter.dump(chr_tx_stat_file)
                    model_reads_printer.flush()
            finally:
                tmp_gff.close()
                tmp_ext_gff.close()
                model_reads_printer.close()
                sqanti_t2t_printer.close()

        open(lock_file, "w", encoding="utf-8").close()
        return dict(aggregator.read_stat_counter.stats_dict), dict(tx_stat_counter.stats_dict)
    finally:
        if genedb is not None:
            genedb.close()


def _filter_umis_worker(task):
    """Worker-logging-safe UMI deduplication pass per chromosome."""
    sample, chr_id, chr_ids, args, edit_distance, output_filtered_reads, barcode_remap, out_prefix_override = task
    out_prefix = out_prefix_override if out_prefix_override else sample.out_umi_filtered_done
    umi_done = umi_filtered_lock_file_name(out_prefix, chr_id, edit_distance)
    allinfo_file = allinfo_file_name(out_prefix, chr_id, edit_distance)
    stats_file = allinfo_stats_file_name(out_prefix, chr_id, edit_distance)

    if args.resume and os.path.exists(umi_done):
        return allinfo_file, stats_file, umi_done
    if os.path.exists(umi_done):
        os.remove(umi_done)

    logger.info("Filtering PCR duplicates for chromosome %s", chr_id)
    transcript_type_dict = create_transcript_info_dict(args.genedb, [chr_id]) if args.genedb else {}
    barcode_feature_table: Dict[str, List[str]] = {}
    if getattr(args, "barcode2spot", None):
        fname, bc_col, spot_cols = parse_barcode2spot_spec(args.barcode2spot)
        with open(fname, "r", encoding="utf-8") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split("\t")
                if len(parts) > bc_col:
                    barcode_feature_table[parts[bc_col]] = [parts[c] for c in spot_cols if c < len(parts)]

    genedb = MmapAnnotationIndex(args.genedb) if args.genedb else None
    try:
        string_pools = setup_string_pools(
            args, sample, chr_ids, chr_id, genedb_index=genedb, load_barcode_pool=True, load_tsv_pools=False
        )
        _load_dynamic_pools(string_pools, sample.get_dynamic_pools_file(chr_id))
        if barcode_remap:
            output_filtered_reads = False
        umi_filter = UMIFilter(args.umi_length, edit_distance, barcode_remap=barcode_remap)
        filtered_reads = sample.get_filtered_reads_file(chr_id) if output_filtered_reads else None
        umi_filter.process_single_chr(
            chr_id, sample.out_raw_file, transcript_type_dict, barcode_feature_table,
            allinfo_file, filtered_reads, stats_file, string_pools,
        )
        open(umi_done, "w", encoding="utf-8").close()
        logger.info("PCR duplicates filtered for chromosome %s", chr_id)
        return allinfo_file, stats_file, umi_done
    finally:
        if genedb is not None:
            genedb.close()


class DatasetProcessor:
    """Orchestrates single-pass O(Locus) streaming and two-stage discovery/SC/multimap processing."""

    def __init__(self, args) -> None:
        self.args = args
        cmd_line = getattr(args, "_cmd_line", "")
        version = getattr(args, "_version", "4.0.0")
        self.common_header = f"# Command line: {cmd_line}\n# IsoQuant version: {version}\n"
        self.grouping_strategy_names = get_grouping_strategy_names(args)
        self.all_read_groups: List[set] = []
        self.alignment_stat_counter = EnumStats()
        self.chr_ids: List[str] = []

        ref_dir = os.path.dirname(args.fai_file_name) if getattr(args, "fai_file_name", None) else None
        with ReferenceGenome(args.reference, index_dir=ref_dir) as ref:
            self.contig_lengths = ref.get_contig_lengths()

    def process_all_samples(self, input_data) -> None:
        logger.info("Running IsoQuant version %s", getattr(self.args, "_version", "4.0.0"))
        for sample in input_data.samples:
            self.process_sample(sample)

    def _can_run_single_pass(self) -> bool:
        return (
            not self.args.use_secondary
            and not self.args.mode.needs_pcr_deduplication()
            and not getattr(self.args, "read_assignments", None)
            and (self.args.no_model_construction or self.args.polya_requirement_strategy != PolyAUsageStrategies.auto)
            and not getattr(self.args, "keep_tmp", False)
        )

    def keep_only_defined_chromosomes(self, chr_set: set) -> set:
        if getattr(self.args, "process_only_chr", None):
            chr_set.intersection_update(self.args.process_only_chr)
        elif getattr(self.args, "discard_chr", None):
            chr_set.difference_update(self.args.discard_chr)
        return chr_set

    def get_chromosome_ids(self, sample: SampleData) -> List[str]:
        genome_chroms = self.keep_only_defined_chromosomes(set(self.contig_lengths.keys()))
        bam_chroms = set()
        for lib in sample.file_list:
            with pysam.AlignmentFile(lib[0], "rb", require_index=True) as af:
                bam_chroms.update(af.references)
        bam_chroms = self.keep_only_defined_chromosomes(bam_chroms)

        overlap = genome_chroms & bam_chroms
        if not overlap:
            raise IsoQuantError(
                IsoQuantExitCode.CHROMOSOME_MISMATCH,
                "Chromosomes in BAM(s) do not overlap reference genome",
            )
        if not self.args.genedb:
            return sorted(overlap, key=lambda c: self.contig_lengths.get(c, 0), reverse=True)

        with MmapAnnotationIndex(self.args.genedb) as idx:
            ann_chroms = set()
            for chrom in idx.seqids():
                g0, g1 = idx.chrom_gene_slice(chrom)
                if g0 < g1:
                    ann_chroms.add(chrom)
        ann_chroms = self.keep_only_defined_chromosomes(ann_chroms)
        common = overlap & ann_chroms
        if not common:
            raise IsoQuantError(
                IsoQuantExitCode.CHROMOSOME_MISMATCH,
                "Chromosomes in annotation do not overlap BAM/reference",
            )
        return sorted(common, key=lambda c: self.contig_lengths.get(c, 0), reverse=True)

    def _prepare_sample_tables(self, sample: SampleData) -> tuple[dict, Optional[str]]:
        self.all_read_groups = [set() for _ in self.grouping_strategy_names]
        rg_lock = read_group_lock_filename(sample)
        if self.args.resume and os.path.exists(rg_lock):
            logger.info("Read group table already split; reusing")
        else:
            if os.path.exists(rg_lock):
                os.remove(rg_lock)
            prepare_read_groups(self.args, sample)
            open(rg_lock, "w", encoding="utf-8").close()

        split_barcodes_dict: Dict[str, str] = {}
        bc_split_done: Optional[str] = None
        if self.args.mode.needs_pcr_deduplication() and not getattr(self.args, "barcoded_bam", False):
            if self.args.barcoded_reads:
                sample.barcoded_reads = list(self.args.barcoded_reads)
            for c in self.chr_ids:
                split_barcodes_dict[c] = f"{sample.barcodes_split_reads}_{convert_chr_id_to_file_name_str(c)}"
            bc_split_done = split_barcodes_lock_filename(sample)
            if not (self.args.resume and os.path.exists(bc_split_done)):
                if os.path.exists(bc_split_done):
                    os.remove(bc_split_done)
                self._split_read_barcode_table(sample)
                open(bc_split_done, "w", encoding="utf-8").close()
        return split_barcodes_dict, bc_split_done

    def _split_read_barcode_table(self, sample: SampleData) -> None:
        if not sample.barcoded_reads:
            return
        bam_file = sample.file_list[0][0]
        # Merge multiple barcode TSVs cleanly if needed before splitting
        for bc_file in sample.barcoded_reads:
            out_map = split_read_table_by_chromosomes(
                bc_file, bam_file, sample.barcodes_split_reads, self.chr_ids, has_header=True
            )
            for chr_id, written in out_map.items():
                expected = sample.get_barcodes_split_file(chr_id)
                if written != expected and os.path.exists(written):
                    if os.path.exists(expected):
                        with open(written, "r", encoding="utf-8") as src, open(expected, "a", encoding="utf-8") as dst:
                            shutil.copyfileobj(src, dst)
                        os.remove(written)
                    else:
                        os.replace(written, expected)

    def _cleanup_split_barcodes(self, split_barcodes_dict: Dict[str, str], bc_split_done: Optional[str]) -> None:
        for path in split_barcodes_dict.values():
            if os.path.exists(path):
                os.remove(path)
        if bc_split_done and os.path.exists(bc_split_done):
            os.remove(bc_split_done)

    def _add_unmapped_counts(self, sample: SampleData) -> None:
        for lib in sample.file_list:
            with pysam.AlignmentFile(lib[0], "rb", require_index=True) as af:
                self.alignment_stat_counter.add(AlignmentType.unaligned, af.unmapped)

    def process_sample(self, sample: SampleData) -> None:
        os.makedirs(sample.out_dir, exist_ok=True)
        os.makedirs(sample.aux_dir, exist_ok=True)
        logger.info("Processing experiment %s", sample.prefix)
        self.chr_ids = self.get_chromosome_ids(sample)
        self.alignment_stat_counter = EnumStats()
        sample.use_technical_replicas = (
            self.args.use_replicas
            and len(sample.file_list) > 1
            and self.args.read_group is not None
            and "file_name" in self.args.read_group
        )

        split_bc_dict, bc_split_done = self._prepare_sample_tables(sample)

        if self._can_run_single_pass():
            self._process_sample_single_pass(sample)
        else:
            if self.args.read_assignments:
                saves_file = self.args.read_assignments[0]
            else:
                self.collect_reads(sample)
                saves_file = sample.out_raw_file
                if self.args.mode.needs_pcr_deduplication():
                    self.filter_umis(sample)

            total_assignments, polya_found, self.all_read_groups = self.load_read_info(saves_file)
            polya_fraction = polya_found / total_assignments if total_assignments > 0 else 0.0
            self._resolve_polya_construction_flags(polya_fraction)
            self.process_assigned_reads(sample, saves_file)

        self._cleanup_split_barcodes(split_bc_dict, bc_split_done)
        logger.info("Processed experiment %s", sample.prefix)

    def _resolve_polya_construction_flags(self, polya_fraction: float) -> None:
        self.args.requires_polya_for_construction = set_polya_requirement_strategy(
            polya_fraction >= self.args.polya_percentage_threshold, self.args.polya_requirement_strategy
        )
        self.args.require_monointronic_polya = set_polya_requirement_strategy(
            self.args.require_monointronic_polya or self.args.requires_polya_for_construction,
            self.args.polya_requirement_strategy,
        )
        self.args.require_monoexonic_polya = set_polya_requirement_strategy(
            self.args.require_monoexonic_polya or self.args.requires_polya_for_construction,
            self.args.polya_requirement_strategy,
        )

    def _process_sample_single_pass(self, sample: SampleData) -> None:
        logger.info("Using single-pass O(Locus) streaming closure for %s", sample.prefix)
        self._resolve_polya_construction_flags(0.0)
        tasks = [(sample, chr_id, self.chr_ids, self.args, self.all_read_groups) for chr_id in self.chr_ids]
        gc.collect()
        with WorkerPoolExecutor(
            max_workers=self.args.threads, debug=getattr(self.args, "debug", False), log_file=getattr(self.args, "log_file", None)
        ) as pool:
            results = pool.map_ordered(_single_pass_chr_worker, tasks)

        with MmapAnnotationIndex(self.args.genedb) if self.args.genedb else _null_cm() as genedb:
            sp = setup_string_pools(self.args, sample, self.chr_ids, genedb_index=genedb)
            with ReadAssignmentAggregator(
                self.args, sample, sp, gzipped=self.args.gzipped, grouping_strategy_names=self.grouping_strategy_names
            ) as aggregator:
                tx_stat_counter = EnumStats()
                for _chr_id, rg, al_stats, read_stats, tx_stats in results:
                    if isinstance(rg, list):
                        for i, gs in enumerate(rg):
                            self.all_read_groups[i].update(gs)
                    elif rg:
                        self.all_read_groups[0].update(rg)
                    for k, v in al_stats.items():
                        self.alignment_stat_counter.add(k, v)
                    for k, v in read_stats.items():
                        aggregator.read_stat_counter.add(k, v)
                    if not self.args.no_model_construction:
                        for k, v in tx_stats.items():
                            tx_stat_counter.add(k, v)

                self._add_unmapped_counts(sample)
                self._merge_outputs_and_finalize(sample, aggregator, tx_stat_counter, self.chr_ids)

    def collect_reads(self, sample: SampleData) -> None:
        logger.info("Collecting read alignments")
        chr_ids = self.chr_ids
        info_file = sample.get_info_file()
        lock_file = sample.get_collection_lock_file()
        if os.path.exists(lock_file):
            if self.args.resume:
                logger.info("Collected reads detected, skipping collection")
                return
            os.remove(lock_file)

        tasks = [(sample, chr_id, chr_ids, self.args) for chr_id in chr_ids]
        gc.collect()
        with WorkerPoolExecutor(
            max_workers=self.args.threads, debug=getattr(self.args, "debug", False), log_file=getattr(self.args, "log_file", None)
        ) as pool:
            results = pool.map_ordered(_collect_reads_worker, tasks)

        total_assignments = 0
        polya_assignments = 0
        need_cross_chr_mm = bool(self.args.use_secondary)
        compact_by_read: Dict[str, List[CompactMultimapTuple]] = defaultdict(list)
        read_counts: Dict[str, int] = defaultdict(int)

        for chr_id, rg, al_stats, _resumed, payloads, cnt_items, l_tot, l_pa in results:
            if isinstance(rg, list):
                for i, gs in enumerate(rg):
                    self.all_read_groups[i].update(gs)
            else:
                self.all_read_groups[0].update(rg)
            for k, v in al_stats.items():
                self.alignment_stat_counter.add(k, v)
            total_assignments += l_tot
            polya_assignments += l_pa
            for b in payloads:
                c = CompactMultimapTuple.from_bytes(b)
                compact_by_read[c.read_id].append(c)
            for rid, cnt in cnt_items:
                read_counts[rid] += cnt

        if need_cross_chr_mm:
            with MmapAnnotationIndex(self.args.genedb) if self.args.genedb else _null_cm() as genedb:
                sp = setup_string_pools(self.args, sample, chr_ids, genedb_index=genedb)
                if not getattr(self.args, "high_memory", False):
                    mm_dict, u_tot, u_pa = prepare_compact_multimappers_from_saves(chr_ids, sample, read_counts, sp)
                else:
                    mm_dict, u_tot, u_pa = compact_by_read, 0, 0
                mm_tot, mm_pa = resolve_compact_multimapper_records(
                    chr_ids, sample, mm_dict, self.args.multimap_strategy, sp
                )
                total_assignments = u_tot + mm_tot
                polya_assignments = u_pa + mm_pa

        self._add_unmapped_counts(sample)
        tmp_info = info_file + f".tmp.{os.getpid()}"
        with open(tmp_info, "wb") as fh:
            write_int(total_assignments, fh)
            write_int(polya_assignments, fh)
            write_int(len(self.all_read_groups), fh)
            for gs in self.all_read_groups:
                write_list(list(gs), fh, write_string)
        os.replace(tmp_info, info_file)
        open(lock_file, "w", encoding="utf-8").close()

    def process_assigned_reads(self, sample: SampleData, dump_filename: str) -> None:
        chr_ids = self.chr_ids
        with MmapAnnotationIndex(self.args.genedb) if self.args.genedb else _null_cm() as genedb:
            sp = (
                setup_string_pools(
                    self.args, sample, chr_ids, chr_id=None, genedb_index=genedb, load_barcode_pool=False, load_tsv_pools=False
                )
                if self.args.read_group
                else None
            )
            with ReadAssignmentAggregator(
                self.args, sample, sp, gzipped=self.args.gzipped, grouping_strategy_names=self.grouping_strategy_names
            ) as aggregator:
                tasks = [(sample, chr_id, chr_ids, dump_filename, self.args, self.all_read_groups) for chr_id in chr_ids]
                gc.collect()
                with WorkerPoolExecutor(
                    max_workers=self.args.threads, debug=getattr(self.args, "debug", False), log_file=getattr(self.args, "log_file", None)
                ) as pool:
                    results = pool.map_ordered(_construct_models_worker, tasks)

                tx_stat_counter = EnumStats()
                for read_stats, tx_stats in results:
                    for k, v in read_stats.items():
                        aggregator.read_stat_counter.add(k, v)
                    if not self.args.no_model_construction:
                        for k, v in tx_stats.items():
                            tx_stat_counter.add(k, v)

                self._merge_outputs_and_finalize(sample, aggregator, tx_stat_counter, chr_ids)

    def _merge_outputs_and_finalize(
        self, sample: SampleData, aggregator: ReadAssignmentAggregator, tx_stat_counter: EnumStats, chr_ids: List[str]
    ) -> None:
        unaligned = self.alignment_stat_counter.stats_dict.get(AlignmentType.unaligned, 0)
        if self.args.genedb and aggregator.read_info_printer:
            merge_files(sample.out_read_info_tsv, sample.prefix, chr_ids, aggregator.read_info_printer.output_file, copy_header=False, header_lines=3)
        if self.args.genedb and aggregator.basic_printer:
            merge_files(sample.out_assigned_tsv, sample.prefix, chr_ids, aggregator.basic_printer.output_file, copy_header=False, header_lines=3)
        if aggregator.corrected_bed_printer:
            merge_files(sample.out_corrected_bed, sample.prefix, chr_ids, aggregator.corrected_bed_printer.output_file, copy_header=False)

        for counter in aggregator.global_counter.counters:
            merge_counts(counter, sample.prefix, chr_ids, unaligned)

        polya_csv = getattr(self.args, "collect_polya_training", None)
        if polya_csv:
            with open(polya_csv, "w", encoding="utf-8") as fh:
                merge_files(f"{sample.out_polya_prediction_tsv}.training.csv", sample.prefix, chr_ids, fh, header_lines=1)
        tss_csv = getattr(self.args, "collect_tss_training", None)
        if tss_csv and self.args.fl_data:
            with open(tss_csv, "w", encoding="utf-8") as fh:
                merge_files(f"{sample.out_tss_prediction_tsv}.training.csv", sample.prefix, chr_ids, fh, header_lines=1)

        if self.args.sqanti_output:
            merge_files(sample.out_t2t_tsv, sample.prefix, chr_ids, aggregator.t2t_sqanti_printer.output_file, copy_header=False, header_lines=1)

        aggregator.global_counter.finalize(self.args)

        if not self.args.no_model_construction:
            exon_id_storage = FeatureIdStorage(SimpleIDDistributor())
            gff_printer = GFFPrinter(sample.out_dir, sample.prefix, exon_id_storage, header=self.common_header)
            model_reads_printer = (
                ReadInfoPrinter(
                    sample.out_transcript_model_reads_tsv, self.args, IOSupport(self.args),
                    additional_header=self.common_header, gzipped=self.args.gzipped,
                )
                if large_output_enabled(self.args, "read2transcripts")
                else VoidPrinter()
            )
            try:
                merge_files(gff_printer.model_fname, sample.prefix, chr_ids, gff_printer.out_gff, copy_header=False)
                if large_output_enabled(self.args, "read2transcripts"):
                    merge_files(
                        sample.out_transcript_model_reads_tsv, sample.prefix, chr_ids,
                        model_reads_printer.output_file, copy_header=False, header_lines=3,
                    )
                for counter in aggregator.transcript_model_global_counter.counters:
                    merge_counts(counter, sample.prefix, chr_ids, unaligned)
                for counter in aggregator.gene_model_global_counter.counters:
                    merge_counts(counter, sample.prefix, chr_ids)
                if self.args.genedb:
                    ext_gff = GFFPrinter(
                        sample.out_dir, sample.prefix, exon_id_storage,
                        gtf_suffix=".extended_annotation.gtf", header=self.common_header,
                    )
                    try:
                        merge_files(ext_gff.model_fname, sample.prefix, chr_ids, ext_gff.out_gff, copy_header=False)
                    finally:
                        ext_gff.close()
                aggregator.transcript_model_global_counter.finalize(self.args)
                aggregator.gene_model_global_counter.finalize(self.args)
            finally:
                gff_printer.close()
                model_reads_printer.close()

        self._cleanup_per_chr_temp_files(sample, chr_ids)

    def _cleanup_per_chr_temp_files(self, sample: SampleData, chr_ids: List[str]) -> None:
        if getattr(self.args, "keep_tmp", False):
            return
        for chr_id in chr_ids:
            chr_prefix = sample.get_chr_prefix(chr_id)
            for f in glob.glob(os.path.join(sample.out_dir, f"{chr_prefix}.*")):
                os.remove(f)

    def filter_umis(self, sample: SampleData) -> None:
        umi_done = umi_filtered_global_lock_file_name(sample.out_umi_filtered_done)
        if self.args.resume and os.path.exists(umi_done):
            return
        if os.path.exists(umi_done):
            os.remove(umi_done)

        for i, ed in enumerate(_UMI_ED_BY_MODE.get(self.args.mode, [])):
            umi_ed_done = umi_filtered_lock_file_name(sample.out_umi_filtered_done, "", ed)
            if self.args.resume and os.path.exists(umi_ed_done):
                return
            if os.path.exists(umi_ed_done):
                os.remove(umi_ed_done)
            out_prefix = umi_output_prefix(sample.out_umi_filtered, ed)
            self._run_umi_round(sample, ed, i == 0, out_prefix, None, None)
            open(umi_ed_done, "w", encoding="utf-8").close()

        if getattr(self.args, "barcode2barcode", None):
            filename, bc_col, spot_cols = parse_barcode2spot_spec(self.args.barcode2barcode)
            full_map = load_barcode2barcode_mapping(filename, bc_col, spot_cols)
            for col_idx in range(len(spot_cols)):
                remap = {bc: spots[col_idx] for bc, spots in full_map.items() if col_idx < len(spots) and spots[col_idx]}
                bc2bc_lock = umi_barcode2barcode_global_lock(sample.out_umi_filtered_done, col_idx)
                if self.args.resume and os.path.exists(bc2bc_lock):
                    continue
                if os.path.exists(bc2bc_lock):
                    os.remove(bc2bc_lock)
                bc2bc_pref = umi_barcode2barcode_prefix(sample.out_umi_filtered_done, col_idx)
                for ed in _UMI_ED_BY_MODE.get(self.args.mode, []):
                    out_prefix = f"{umi_output_prefix(sample.out_umi_filtered, ed)}_col{col_idx}"
                    self._run_umi_round(sample, ed, False, out_prefix, remap, bc2bc_pref)
                open(bc2bc_lock, "w", encoding="utf-8").close()

        open(umi_done, "w", encoding="utf-8").close()

    def _run_umi_round(self, sample: SampleData, ed: int, output_filtered: bool, out_prefix: str, remap, out_override) -> None:
        tasks = [
            (sample, chr_id, self.chr_ids, self.args, ed, output_filtered, remap, out_override)
            for chr_id in self.chr_ids
        ]
        with WorkerPoolExecutor(
            max_workers=self.args.threads, debug=getattr(self.args, "debug", False), log_file=getattr(self.args, "log_file", None)
        ) as pool:
            results = pool.map_ordered(_filter_umis_worker, tasks)

        save_allinfo = large_output_enabled(self.args, "allinfo")
        allinfo_fh = None
        if save_allinfo:
            fname = f"{out_prefix}.allinfo" + (".gz" if self.args.gzipped else "")
            allinfo_fh = gzip.open(fname, "wt") if self.args.gzipped else open(fname, "w", encoding="utf-8")

        stat_dict: Dict[str, int] = defaultdict(int)
        to_remove: List[str] = []
        try:
            for all_f, st_f, done_f in results:
                if allinfo_fh is not None and os.path.exists(all_f):
                    with open(all_f, "r", encoding="utf-8") as src:
                        shutil.copyfileobj(src, allinfo_fh)
                if os.path.exists(st_f):
                    with open(st_f, "r", encoding="utf-8") as src:
                        for ln in src:
                            parts = ln.strip().split("\t")
                            if len(parts) == 2:
                                stat_dict[parts[0]] += int(parts[1])
                to_remove.extend([all_f, st_f, done_f])
        finally:
            if allinfo_fh is not None:
                allinfo_fh.close()

        with open(f"{out_prefix}.stats.tsv", "w", encoding="utf-8") as outf:
            for k, v in stat_dict.items():
                outf.write(f"{k}\t{v}\n")
        for f in to_remove:
            if os.path.exists(f):
                os.remove(f)

    @staticmethod
    def load_read_info(dump_filename: str):
        with open(info_file_name(dump_filename), "rb") as fh:
            total = read_int(fh)
            polya = read_int(fh)
            n_strat = read_int(fh)
            groups = [set(read_list(fh, read_string)) for _ in range(n_strat)]
        return total, polya, groups


class _null_cm:
    def __enter__(self):
        return None

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        return None
