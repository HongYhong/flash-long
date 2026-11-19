"""UMI deduplication engine using Cython Levenshtein distance with length-difference pruning."""

from __future__ import annotations

import logging
import os
import sqlite3
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import pysam

from isoquant_next.barcodes.common import edit_distance_fast

logger = logging.getLogger("IsoQuant")

UNTRUSTED_UMI_KEY = None
UNTRUSTED_UMIS = frozenset({"", "None", None, "*", "."})


def is_untrusted_umi(umi: Optional[str]) -> bool:
    return umi in UNTRUSTED_UMIS


def format_read_assignment_for_output(read_assignment) -> str:
    """Format ReadAssignment record into .allinfo TSV line safely without attribute mutations."""
    gene_id = "None"
    transcript_id = "None"
    matches = getattr(read_assignment, "isoform_matches", None)
    if matches:
        transcript_id = matches[0].assigned_transcript or "None"
        gene_id = matches[0].assigned_gene or "None"
    attrs = getattr(read_assignment, "additional_attributes", {}) or {}
    t_type = attrs.get("transcript_type", "unknown")
    polya_site = attrs.get("polya_site", -1)
    strand = getattr(read_assignment, "strand", ".")
    barcode = getattr(read_assignment, "barcode", None) or "*"
    umi = getattr(read_assignment, "umi", None) or "*"
    exons_str = ",".join(f"{s}-{e}" for s, e in getattr(read_assignment, "corrected_exons", []))
    return (
        f"{read_assignment.read_id}\t{gene_id}\t{transcript_id}\t"
        f"{read_assignment.chr_id}:{read_assignment.start}-{read_assignment.end}\t"
        f"{strand}\t{barcode}\t{umi}\t{t_type}\t{polya_site}\t{exons_str}"
    )


class UMIFilter:
    """Fast length-pruned Cython Levenshtein UMI cluster deduplicator."""

    def __init__(
        self,
        umi_length: int = 0,
        edit_distance: int = 3,
        disregard_length_diff: bool = True,
        only_unique_assignments: bool = False,
        only_spliced_reads: bool = False,
        barcode_remap: Optional[Dict[str, str]] = None,
    ) -> None:
        self.umi_length = int(umi_length)
        self.max_edit_distance = int(edit_distance)
        self.disregard_length_diff = bool(disregard_length_diff)
        self.only_unique_assignments = bool(only_unique_assignments)
        self.only_spliced_reads = bool(only_spliced_reads)
        self.barcode_remap = barcode_remap or {}
        self.stats: Dict[str, int] = defaultdict(int)
        self.duplicated_molecule_counts: Dict[int, int] = defaultdict(int)
        self.unique_gene_barcode: Set = set()
        self.selected_reads: Set[str] = set()
        self.total_assignments: int = 0

    def _find_similar_umi(self, umi: str, trusted_umi_list: List[str]) -> Optional[str]:
        if self.max_edit_distance == -1:
            return None if not trusted_umi_list else trusted_umi_list[0]

        similar_umi: Optional[str] = None
        best_dist = 1 << 30
        max_ed = self.max_edit_distance
        disregard = self.disregard_length_diff
        lu = len(umi)

        for occ in trusted_umi_list:
            if max_ed == 0:
                if disregard:
                    similar, ed = occ == umi, 0
                elif len(occ) < lu:
                    similar, ed = occ in umi, lu - len(occ)
                else:
                    similar, ed = umi in occ, len(occ) - lu
            elif occ == umi:
                similar, ed = True, 0
            else:
                len_diff = abs(len(occ) - lu)
                if disregard and len_diff > max_ed:
                    # O(1) length-difference pruning avoids entering C DP when unadjusted ed cannot satisfy <= max_ed
                    continue
                ed = edit_distance_fast(occ, umi)
                if not disregard:
                    ed -= len_diff
                similar = ed <= max_ed

            if similar and ed < best_dist:
                similar_umi = occ
                best_dist = ed
                if best_dist == 0:
                    break
        return similar_umi

    def _construct_umi_dict(self, molecule_list: List) -> Dict[Optional[str], List]:
        trusted: Dict[str, List] = {}
        untrusted: List = []
        # Sort by UMI length closeness to expected length, then lexicographically for determinism
        sorted_mols = sorted(
            molecule_list,
            key=lambda m: (is_untrusted_umi(m.umi), abs(len(m.umi or "") - self.umi_length), m.umi or "", m.read_id),
        )
        for m in sorted_mols:
            if is_untrusted_umi(m.umi):
                untrusted.append(m)
                continue
            umi = m.umi
            assert umi is not None
            match = self._find_similar_umi(umi, list(trusted.keys()))
            if match is None:
                trusted[umi] = [m]
            else:
                trusted[match].append(m)

        out: Dict[Optional[str], List] = dict(trusted)
        if untrusted:
            out[UNTRUSTED_UMI_KEY] = untrusted
        return out

    def _process_duplicates(self, molecule_list: List) -> List:
        if not molecule_list:
            return []
        if len(molecule_list) == 1:
            self.duplicated_molecule_counts[1] += 1
            return list(molecule_list)

        resulting_reads = []
        umi_dict = self._construct_umi_dict(molecule_list)

        for umi, dups in umi_dict.items():
            self.duplicated_molecule_counts[len(dups)] += 1
            if len(dups) == 1:
                resulting_reads.append((dups[0], umi))
                continue

            best_read = dups[0]
            for m in dups[1:]:
                if not best_read.assignment_type.is_unique() and m.assignment_type.is_unique():
                    best_read = m
                elif len(m.corrected_exons) > len(best_read.corrected_exons):
                    best_read = m
                elif (
                    len(m.corrected_exons) == len(best_read.corrected_exons)
                    and (m.corrected_exons[-1][1] - m.corrected_exons[0][0])
                    > (best_read.corrected_exons[-1][1] - best_read.corrected_exons[0][0])
                ):
                    best_read = m

            polyas = set()
            transcript_types = set()
            isoform_ids = set()
            for m in dups:
                if m.read_id != best_read.read_id:
                    continue
                transcript_types.add(m.additional_attributes.get("transcript_type", "unknown"))
                polyas.add(m.additional_attributes.get("polya_site", -1))
                if m.isoform_matches:
                    isoform_ids.add(m.isoform_matches[0].assigned_transcript)

            if len(transcript_types) > 1:
                best_read.set_additional_attribute("transcript_type", "None")
            if len(polyas) > 1:
                best_read.set_additional_attribute("polya_site", -1)
            if len(isoform_ids) > 1 and best_read.isoform_matches:
                best_read.isoform_matches[0].assigned_transcript = "None"

            if best_read.assignment_type.is_inconsistent():
                if best_read.isoform_matches:
                    best_read.isoform_matches[0].assigned_transcript = "None"
                best_read.set_additional_attribute("polya_site", -1)
                best_read.set_additional_attribute("transcript_type", "None")

            resulting_reads.append((best_read, umi))

        if len(resulting_reads) == 1:
            return [resulting_reads[0][0]]
        return [x[0] for x in resulting_reads if x[1] != UNTRUSTED_UMI_KEY]

    def _process_gene(self, gene_dict: Dict[str, List]) -> List:
        kept: List = []
        for _barcode, reads in gene_dict.items():
            if not reads:
                continue
            kept.extend(self._process_duplicates(reads))
        return kept

    def _process_chunk(
        self,
        gene_barcode_dict: Dict[str, Dict[str, List]],
        allinfo_outf,
        read_ids_outf=None,
    ) -> Tuple[int, int]:
        total_in = 0
        total_kept = 0
        for gene_id, bc_map in gene_barcode_dict.items():
            for _bc, reads in bc_map.items():
                total_in += len(reads)
            kept = self._process_gene(bc_map)
            total_kept += len(kept)
            for r in kept:
                self.add_stats_for_read(r)
                self.selected_reads.add(r.read_id)
                if allinfo_outf is not None:
                    allinfo_outf.write(format_read_assignment_for_output(r) + "\n")
                if read_ids_outf is not None:
                    read_ids_outf.write(r.read_id + "\n")
        return total_in, total_kept

    def add_stats_for_read(self, read_assignment) -> None:
        self.stats["selected_reads"] += 1
        barcode = getattr(read_assignment, "barcode", None)
        gene_id = None
        matches = getattr(read_assignment, "isoform_matches", None)
        if matches:
            gene_id = matches[0].assigned_gene
        if gene_id and barcode:
            self.unique_gene_barcode.add((gene_id, barcode))

    def process_single_chr(
        self,
        chr_id: str,
        saves_prefix: str,
        transcript_type_dict: Dict[str, Tuple[str, int]],
        barcode_feature_table: Dict[str, str],
        all_info_file_name: str,
        filtered_reads_file_name: Optional[str],
        stats_output_file_name: str,
        string_pools=None,
    ) -> Tuple[str, str]:
        from isoquant_next.assignment.loader import create_merging_assignment_loader

        os.makedirs(os.path.dirname(all_info_file_name) or ".", exist_ok=True)
        read_count = 0
        spliced_count = 0
        self.unique_gene_barcode = set()

        with open(all_info_file_name, "w", encoding="utf-8") as allinfo_outf:
            filtered_reads_outf = open(filtered_reads_file_name, "w", encoding="utf-8") if filtered_reads_file_name else None
            try:
                with create_merging_assignment_loader(chr_id, saves_prefix, string_pools=string_pools) as loader:
                    while loader.has_next():
                        gene_barcode_dict: Dict[str, Dict[str, List]] = defaultdict(lambda: defaultdict(list))
                        _, assignment_storage = loader.get_next()
                        if not assignment_storage:
                            continue
                        for read_assignment in assignment_storage:
                            if read_assignment.gene_assignment_type.is_unassigned():
                                continue
                            if read_assignment.gene_assignment_type.is_ambiguous():
                                self.stats["Assigned to any gene"] += 1
                                if len(read_assignment.corrected_exons) > 1:
                                    self.stats["Spliced"] += 1
                                continue

                            barcode = read_assignment.barcode
                            if barcode == "*" or not barcode:
                                barcode = None
                            spliced = len(read_assignment.corrected_exons) > 1
                            barcoded = barcode is not None

                            if barcode is None or barcode not in barcode_feature_table:
                                cell_type = "None"
                            else:
                                cell_types = barcode_feature_table[barcode]
                                cell_type = ",".join(cell_types) if isinstance(cell_types, list) and cell_types else (
                                    "None" if isinstance(cell_types, list) else cell_types
                                )
                            read_assignment.set_additional_attribute("cell_type", cell_type)

                            if len(read_assignment.isoform_matches) == 1:
                                isoform_match = read_assignment.isoform_matches[0]
                                if read_assignment.assignment_type.is_consistent():
                                    transcript_id = isoform_match.assigned_transcript
                                    t_type, polya_site = transcript_type_dict.get(transcript_id, ("unknown_type", -1))
                                else:
                                    t_type, polya_site = "unknown_type", -1
                                read_assignment.set_additional_attribute("transcript_type", t_type)
                                read_assignment.set_additional_attribute("polya_site", polya_site)
                            elif read_assignment.assignment_type.is_consistent():
                                t_types = set()
                                polya_sites = set()
                                for m in read_assignment.isoform_matches:
                                    t_type, polya_site = transcript_type_dict.get(m.assigned_transcript, ("unknown_type", -1))
                                    t_types.add(t_type)
                                    polya_sites.add(polya_site)
                                read_assignment.set_additional_attribute(
                                    "transcript_type", t_types.pop() if len(t_types) == 1 else "None"
                                )
                                read_assignment.set_additional_attribute(
                                    "polya_site", polya_sites.pop() if len(polya_sites) == 1 else -1
                                )
                                if read_assignment.isoform_matches:
                                    read_assignment.isoform_matches[0].assigned_transcript = "None"
                            else:
                                self.total_assignments += 1
                                continue

                            self.total_assignments += 1
                            self.add_stats_for_read(read_assignment)
                            if not barcoded or (not spliced and self.only_spliced_reads):
                                continue
                            if read_assignment.isoform_matches:
                                gene_id = read_assignment.isoform_matches[0].assigned_gene
                                grouping_key = barcode
                                if self.barcode_remap and barcode in self.barcode_remap:
                                    grouping_key = self.barcode_remap[barcode]
                                    read_assignment.barcode = grouping_key
                                gene_barcode_dict[gene_id][grouping_key].append(read_assignment)

                        chunk_in, chunk_kept = self._process_chunk(gene_barcode_dict, allinfo_outf, filtered_reads_outf)
                        read_count += chunk_in
                        spliced_count += chunk_kept
            finally:
                if filtered_reads_outf is not None:
                    filtered_reads_outf.close()

        with open(stats_output_file_name, "w", encoding="utf-8") as fh:
            fh.write(f"Unique gene-barcodes pairs\t{len(self.unique_gene_barcode)}\n")
            fh.write(f"Total reads saved\t{read_count}\n")
            fh.write(f"Spliced reads saved\t{spliced_count}\n")
            fh.write(f"Total assignments processed\t{self.total_assignments}\n")
            for k in sorted(self.stats.keys()):
                fh.write(f"{k}\t{self.stats[k]}\n")
        return all_info_file_name, stats_output_file_name


def filter_bam(in_file_name: str, out_file_name: str, read_set: Set[str]) -> None:
    """Copy BAM alignments whose query_name belongs to read_set."""
    with pysam.AlignmentFile(in_file_name, "rb") as src:
        with pysam.AlignmentFile(out_file_name, "wb", template=src) as dst:
            for aln in src.fetch(until_eof=True):
                if aln.query_name in read_set:
                    dst.write(aln)


def create_transcript_info_dict(genedb: str, chr_ids: Optional[List[str]] = None) -> Dict[str, Tuple[str, int]]:
    """Build {transcript_id: (biotype, exon_count)} from .iqidx (zero-ORM) or gffutils .db."""
    path_str = str(genedb)
    if path_str.endswith(".iqidx"):
        from isoquant_next.annotation.iqidx import MmapAnnotationIndex

        out: Dict[str, Tuple[str, int]] = {}
        target_chroms = set(chr_ids) if chr_ids else None
        with MmapAnnotationIndex(path_str) as idx:
            g_tx = idx.gene_tx_ptr
            t_id_s = idx.tx_id_str
            t_attr_s = idx.tx_attr_str
            t_ex = idx.tx_exon_ptr
            try:
                for chrom in idx.seqids():
                    if target_chroms is not None and chrom not in target_chroms:
                        continue
                    g0, g1 = idx.chrom_gene_slice(chrom)
                    for gi in range(g0, g1):
                        for ti in range(int(g_tx[gi]), int(g_tx[gi + 1])):
                            tx_id = idx.get_str(int(t_id_s[ti]))
                            attr = idx.get_str(int(t_attr_s[ti]))
                            biotype = "protein_coding" if "protein_coding" in attr else "unknown"
                            n_exons = int(t_ex[ti + 1] - t_ex[ti])
                            out[tx_id] = (biotype, n_exons)
            finally:
                del g_tx, t_id_s, t_attr_s, t_ex
        return out

    out = {}
    conn = sqlite3.connect(f"file:{path_str}?mode=ro", uri=True)
    try:
        cur = conn.cursor()
        sql = "SELECT id, seqid, attributes FROM features WHERE featuretype IN ('transcript', 'mRNA')"
        rows = cur.execute(sql).fetchall()
        target_chroms = set(chr_ids) if chr_ids else None
        tx_meta: Dict[str, str] = {}
        for tid, seqid, attrs in rows:
            if target_chroms is not None and seqid not in target_chroms:
                continue
            tx_meta[tid] = "protein_coding" if (attrs and "protein_coding" in attrs) else "unknown"
        rel_rows = cur.execute("SELECT parent, COUNT(*) FROM relations WHERE level = 1 GROUP BY parent").fetchall()
        rel_cnt = {p: int(c) for p, c in rel_rows}
        for tid, bt in tx_meta.items():
            out[tid] = (bt, rel_cnt.get(tid, 1))
    finally:
        conn.close()
    return out
