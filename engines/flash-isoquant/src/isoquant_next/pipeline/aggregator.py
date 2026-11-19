"""ReadAssignmentAggregator wiring Node 4 printers, Node 5 terminal counters, and Node 7 quant counters."""

from __future__ import annotations

import logging
from typing import Optional, Set

from isoquant_next.annotation.iqidx import MmapAnnotationIndex
from isoquant_next.assignment.io import (
    BasicTSVAssignmentPrinter,
    BEDPrinter,
    IOSupport,
    ReadAssignmentCompositePrinter,
    ReadInfoPrinter,
    SqantiTSVPrinter,
)
from isoquant_next.core.stats import EnumStats
from isoquant_next.core.utils import large_output_enabled
from isoquant_next.discovery.printer import VoidTranscriptPrinter
from isoquant_next.quantification import (
    CompositeCounter,
    ExonCounter,
    ExonSpliceSiteCounter,
    IntronCounter,
    IntronRetentionCounter,
    JointExonCounter,
    create_gene_counter,
    create_transcript_counter,
)
from isoquant_next.terminal.counters import PolyACounter, TSSCounter

logger = logging.getLogger("IsoQuant")


def _chr_features_from_index(genedb, chr_id: str) -> tuple[Set[str], Set[str]]:
    if genedb is None or not chr_id:
        return set(), set()
    if isinstance(genedb, MmapAnnotationIndex):
        gene_set: Set[str] = set()
        tx_set: Set[str] = set()
        g0, g1 = genedb.chrom_gene_slice(chr_id)
        g_id_s = genedb.gene_id_str
        g_tx = genedb.gene_tx_ptr
        t_id_s = genedb.tx_id_str
        for gi in range(g0, g1):
            gene_set.add(genedb.get_str(int(g_id_s[gi])))
            for ti in range(int(g_tx[gi]), int(g_tx[gi + 1])):
                tx_set.add(genedb.get_str(int(t_id_s[ti])))
        return gene_set, tx_set

    gene_set = {g.id for g in genedb.region(seqid=chr_id, start=1, featuretype="gene")}
    tx_set = {t.id for t in genedb.region(seqid=chr_id, start=1, featuretype=("transcript", "mRNA"))}
    return gene_set, tx_set


class ReadAssignmentAggregator:
    """Manages printers and counters (ungrouped + grouped + discovery models) for a sample or chromosome."""

    def __init__(
        self,
        args,
        sample,
        string_pools,
        genedb=None,
        chr_id: Optional[str] = None,
        gzipped: bool = False,
        grouping_strategy_names=None,
    ) -> None:
        self.args = args
        self.string_pools = string_pools
        self.grouping_strategy_names = grouping_strategy_names if grouping_strategy_names else ["default"]
        cmd_line = getattr(args, "_cmd_line", "")
        version = getattr(args, "_version", "4.0.0")
        self.common_header = f"# Command line: {cmd_line}\n# IsoQuant version: {version}\n"
        self.io_support = IOSupport(self.args)

        self.gene_set, self.transcript_set = _chr_features_from_index(genedb, chr_id)
        self.read_stat_counter = EnumStats()
        self.polya_counter: Optional[PolyACounter] = None
        self.tss_counter: Optional[TSSCounter] = None
        self.gene_counter = None
        self.transcript_counter = None
        self.transcript_model_counter = None
        self.gene_model_counter = None

        self._init_printers(sample, chr_id, gzipped)
        self._init_ungrouped_counters(sample, chr_id)
        self._init_grouped_counters(sample, chr_id)
        self._init_grouped_model_counters(sample, chr_id)

    def close(self) -> None:
        for p in (self.corrected_bed_printer, self.read_info_printer, self.basic_printer, self.t2t_sqanti_printer):
            if p is not None and hasattr(p, "close"):
                p.close()

    def __enter__(self) -> "ReadAssignmentAggregator":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def _init_printers(self, sample, chr_id, gzipped: bool) -> None:
        printer_list = []
        self.corrected_bed_printer = self._make_corrected_bed_printer(sample, chr_id, gzipped, printer_list)
        self.read_info_printer = self._make_read_info_printer(sample, chr_id, gzipped, printer_list)
        self.basic_printer = self._make_basic_printer(sample, chr_id, gzipped, printer_list)
        self.t2t_sqanti_printer = self._make_t2t_printer(sample, chr_id)
        self.global_printer = ReadAssignmentCompositePrinter(printer_list)

    def _make_corrected_bed_printer(self, sample, chr_id, gzipped: bool, printer_list):
        if not large_output_enabled(self.args, "corrected_bed"):
            return None
        path = sample.get_corrected_bed_file(chr_id) if chr_id else sample.out_corrected_bed
        printer = BEDPrinter(path, self.args, print_corrected=True, gzipped=gzipped)
        printer_list.append(printer)
        return printer

    def _make_read_info_printer(self, sample, chr_id, gzipped: bool, printer_list):
        if not (self.args.genedb and large_output_enabled(self.args, "read_info")):
            return None
        path = sample.get_read_info_tsv_file(chr_id) if chr_id else sample.out_read_info_tsv
        printer = ReadInfoPrinter(path, self.args, self.io_support, additional_header=self.common_header, gzipped=gzipped)
        printer_list.append(printer)
        return printer

    def _make_basic_printer(self, sample, chr_id, gzipped: bool, printer_list):
        if not (self.args.genedb and large_output_enabled(self.args, "read_assignments")):
            return None
        path = sample.get_assigned_tsv_file(chr_id) if chr_id else sample.out_assigned_tsv
        printer = BasicTSVAssignmentPrinter(
            path, self.args, self.io_support, additional_header=self.common_header, gzipped=gzipped
        )
        sample.out_assigned_tsv_result = printer.output_file_name
        printer_list.append(printer)
        return printer

    def _make_t2t_printer(self, sample, chr_id):
        if not getattr(self.args, "sqanti_output", False):
            return VoidTranscriptPrinter()
        path = sample.get_t2t_tsv_file(chr_id) if chr_id else sample.out_t2t_tsv
        return SqantiTSVPrinter(path, self.args, self.io_support)

    def _init_ungrouped_counters(self, sample, chr_id) -> None:
        self.global_counter = CompositeCounter()
        self._add_ungrouped_quant_counters(sample, chr_id)

        self.transcript_model_global_counter = CompositeCounter()
        self.gene_model_global_counter = CompositeCounter()
        self._add_ungrouped_model_counters(sample, chr_id)

        self._add_ungrouped_exon_counters(sample, chr_id)
        self._add_ungrouped_terminal_counters(sample, chr_id)
        self._add_ungrouped_ir_counter(sample, chr_id)

    def _add_ungrouped_quant_counters(self, sample, chr_id) -> None:
        if not (self.args.genedb and self.args.run_quantification):
            return
        gene_path = sample.get_gene_counts_file(chr_id) if chr_id else sample.out_gene_counts_tsv
        tx_path = sample.get_transcript_counts_file(chr_id) if chr_id else sample.out_transcript_counts_tsv
        self.gene_counter = create_gene_counter(
            gene_path, self.args.gene_quantification, complete_feature_list=self.gene_set
        )
        self.transcript_counter = create_transcript_counter(
            tx_path, self.args.transcript_quantification, complete_feature_list=self.transcript_set
        )
        self.global_counter.add_counters([self.gene_counter, self.transcript_counter])

    def _add_ungrouped_model_counters(self, sample, chr_id) -> None:
        if self.args.no_model_construction:
            return
        tx_path = sample.get_transcript_model_counts_file(chr_id) if chr_id else sample.out_transcript_model_counts_tsv
        gene_path = sample.get_gene_model_counts_file(chr_id) if chr_id else sample.out_gene_model_counts_tsv
        self.transcript_model_counter = create_transcript_counter(tx_path, self.args.transcript_quantification)
        self.gene_model_counter = create_gene_counter(gene_path, self.args.gene_quantification)
        self.transcript_model_global_counter.add_counter(self.transcript_model_counter)
        self.gene_model_global_counter.add_counter(self.gene_model_counter)

    def _add_ungrouped_exon_counters(self, sample, chr_id) -> None:
        if not (self.args.count_exons and self.args.genedb):
            return
        exon_path = sample.get_exon_counts_file(chr_id) if chr_id else sample.out_exon_counts_tsv
        intron_path = sample.get_intron_counts_file(chr_id) if chr_id else sample.out_intron_counts_tsv
        ess_path = (
            sample.get_exon_splice_site_counts_file(chr_id)
            if chr_id
            else sample.out_exon_splice_site_counts_tsv
        )
        self.exon_counter = JointExonCounter(exon_path)
        self.intron_counter = IntronCounter(intron_path)
        self.exon_splice_site_counter = ExonSpliceSiteCounter(
            ess_path, delta=self.args.delta, emit_read_ids=getattr(self.args, "emit_read_ids", False)
        )
        counters = [self.exon_counter, self.intron_counter, self.exon_splice_site_counter]
        if getattr(self.args, "old_exon_count_format", False):
            old_path = sample.get_old_exon_counts_file(chr_id) if chr_id else sample.out_old_exon_counts_tsv
            self.old_exon_counter = ExonCounter(old_path)
            counters.append(self.old_exon_counter)
        self.global_counter.add_counters(counters)

    def _add_ungrouped_terminal_counters(self, sample, chr_id) -> None:
        if not self.args.predict_terminal_sites:
            return
        polya_path = sample.get_polya_prediction_file(chr_id) if chr_id else sample.out_polya_prediction_tsv
        self.polya_counter = PolyACounter(self.args, polya_path)
        self.global_counter.add_counter(self.polya_counter)
        if self.args.fl_data:
            tss_path = sample.get_tss_prediction_file(chr_id) if chr_id else sample.out_tss_prediction_tsv
            self.tss_counter = TSSCounter(self.args, tss_path)
            self.global_counter.add_counter(self.tss_counter)

    def _add_ungrouped_ir_counter(self, sample, chr_id) -> None:
        if not (self.args.count_intron_retentions and self.args.genedb):
            return
        ir_path = sample.get_intron_retention_counts_file(chr_id) if chr_id else sample.out_intron_retention_counts_tsv
        self.intron_retention_counter = IntronRetentionCounter(ir_path)
        self.global_counter.add_counter(self.intron_retention_counter)

    def _init_grouped_counters(self, sample, chr_id) -> None:
        if not (self.args.read_group and self.args.genedb):
            return
        for group_idx, strategy_name in enumerate(self.grouping_strategy_names):
            if self.args.run_quantification:
                self._add_grouped_quant_counters(sample, chr_id, group_idx, strategy_name)
            if self.args.count_exons:
                self._add_grouped_exon_counters(sample, chr_id, group_idx, strategy_name)
            self._add_grouped_terminal_counters(sample, chr_id, group_idx, strategy_name)
            if self.args.count_intron_retentions:
                self._add_grouped_ir_counter(sample, chr_id, group_idx, strategy_name)

    def _add_grouped_quant_counters(self, sample, chr_id, group_idx: int, strategy_name: str) -> None:
        gene_out = (
            sample.get_grouped_counts_file(chr_id, "gene", strategy_name)
            if chr_id
            else f"{sample.out_gene_grouped_counts_tsv}_{strategy_name}"
        )
        tx_out = (
            sample.get_grouped_counts_file(chr_id, "transcript", strategy_name)
            if chr_id
            else f"{sample.out_transcript_grouped_counts_tsv}_{strategy_name}"
        )
        self.global_counter.add_counters([
            create_gene_counter(
                gene_out,
                self.args.gene_quantification,
                complete_feature_list=self.gene_set,
                string_pools=self.string_pools,
                group_index=group_idx,
            ),
            create_transcript_counter(
                tx_out,
                self.args.transcript_quantification,
                complete_feature_list=self.transcript_set,
                string_pools=self.string_pools,
                group_index=group_idx,
            ),
        ])

    def _add_grouped_exon_counters(self, sample, chr_id, group_idx: int, strategy_name: str) -> None:
        exon_out = (
            sample.get_grouped_counts_file(chr_id, "exon", strategy_name)
            if chr_id
            else f"{sample.out_exon_grouped_counts_tsv}_{strategy_name}"
        )
        intron_out = (
            sample.get_grouped_counts_file(chr_id, "splice_junction", strategy_name)
            if chr_id
            else f"{sample.out_intron_grouped_counts_tsv}_{strategy_name}"
        )
        ess_out = (
            sample.get_grouped_counts_file(chr_id, "exon_splice_site", strategy_name)
            if chr_id
            else f"{sample.out_exon_splice_site_grouped_counts_tsv}_{strategy_name}"
        )
        counters = [
            JointExonCounter(exon_out, string_pools=self.string_pools, group_index=group_idx),
            IntronCounter(intron_out, string_pools=self.string_pools, group_index=group_idx),
            ExonSpliceSiteCounter(
                ess_out,
                string_pools=self.string_pools,
                group_index=group_idx,
                delta=self.args.delta,
                emit_read_ids=getattr(self.args, "emit_read_ids", False),
            ),
        ]
        if getattr(self.args, "old_exon_count_format", False):
            old_out = (
                sample.get_grouped_counts_file(chr_id, "old_exon", strategy_name)
                if chr_id
                else f"{sample.out_old_exon_grouped_counts_tsv}_{strategy_name}"
            )
            counters.append(ExonCounter(old_out, string_pools=self.string_pools, group_index=group_idx))
        self.global_counter.add_counters(counters)

    def _add_grouped_terminal_counters(self, sample, chr_id, group_idx: int, strategy_name: str) -> None:
        if self.args.predict_terminal_sites and not getattr(self.args, "collect_polya_training", None):
            out_f = (
                sample.get_grouped_counts_file(chr_id, "polyA_prediction", strategy_name) + ".tsv"
                if chr_id
                else f"{sample.out_polya_prediction_grouped_tsv}_{strategy_name}.tsv"
            )
            self.global_counter.add_counter(
                PolyACounter(self.args, out_f, string_pools=self.string_pools, group_index=group_idx)
            )
        if self.args.predict_terminal_sites and self.args.fl_data and not getattr(self.args, "collect_tss_training", None):
            out_f = (
                sample.get_grouped_counts_file(chr_id, "TSS_prediction", strategy_name) + ".tsv"
                if chr_id
                else f"{sample.out_tss_prediction_grouped_tsv}_{strategy_name}.tsv"
            )
            self.global_counter.add_counter(
                TSSCounter(self.args, out_f, string_pools=self.string_pools, group_index=group_idx)
            )

    def _add_grouped_ir_counter(self, sample, chr_id, group_idx: int, strategy_name: str) -> None:
        ir_out = (
            sample.get_grouped_counts_file(chr_id, "intron_retention", strategy_name)
            if chr_id
            else f"{sample.out_intron_retention_grouped_counts_tsv}_{strategy_name}"
        )
        self.global_counter.add_counter(IntronRetentionCounter(ir_out, string_pools=self.string_pools, group_index=group_idx))

    def _init_grouped_model_counters(self, sample, chr_id) -> None:
        if not (self.args.read_group and not self.args.no_model_construction):
            return
        for group_idx, strategy_name in enumerate(self.grouping_strategy_names):
            tx_out = (
                sample.get_grouped_counts_file(chr_id, "discovered_transcript", strategy_name)
                if chr_id
                else f"{sample.out_transcript_model_grouped_counts_tsv}_{strategy_name}"
            )
            gene_out = (
                sample.get_grouped_counts_file(chr_id, "discovered_gene", strategy_name)
                if chr_id
                else f"{sample.out_gene_model_grouped_counts_tsv}_{strategy_name}"
            )
            self.transcript_model_global_counter.add_counter(
                create_transcript_counter(
                    tx_out, self.args.transcript_quantification, string_pools=self.string_pools, group_index=group_idx
                )
            )
            self.gene_model_global_counter.add_counter(
                create_gene_counter(
                    gene_out, self.args.gene_quantification, string_pools=self.string_pools, group_index=group_idx
                )
            )
