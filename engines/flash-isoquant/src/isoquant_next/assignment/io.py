"""Assignment output printers, binary temporary writers/loaders (context-managed), and IOSupport."""

from __future__ import annotations

import gzip
import logging
from typing import Iterable, List, Optional, Sequence, Tuple

from isoquant_next.annotation.gene_info import GeneInfo, GeneList
from isoquant_next.assignment.types import (
    BasicReadAssignment,
    MatchClassification,
    ReadAssignment,
    ReadAssignmentType,
    match_subtype_to_str_with_additional_info,
)
from isoquant_next.core.serialization import SHORT_TERMINATION_INT, read_short_int, write_short_int
from isoquant_next.core.utils import (
    CANONICAL_FWD_SITES,
    CANONICAL_REV_SITES,
    argmin,
    junctions_from_blocks,
    range_list_to_str,
    sum_intervals_from_point,
    sum_intervals_to_point,
)

logger = logging.getLogger("IsoQuant")


class PrintAllFunctor:
    def check(self, assignment) -> bool:
        return True


class PrintOnlyFunctor:
    def __init__(self, allowed_types) -> None:
        if isinstance(allowed_types, (list, set, tuple)):
            self.allowed_types = set(allowed_types)
        else:
            self.allowed_types = {allowed_types}

    def check(self, assignment) -> bool:
        if assignment is None:
            return False
        return assignment.assignment_type in self.allowed_types


class AbstractAssignmentPrinter:
    def __init__(self, output_file_name: Optional[str], params, assignment_checker=None) -> None:
        self.assignment_checker = PrintAllFunctor() if assignment_checker is None else assignment_checker
        self.params = params
        self.output_file_name = output_file_name

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def add_read_info(self, read_assignment) -> None:
        raise NotImplementedError()

    def flush(self) -> None:
        raise NotImplementedError()

    def close(self) -> None:
        pass


class TextFileAssignmentPrinter(AbstractAssignmentPrinter):
    def __init__(
        self,
        output_file_name: str,
        params,
        assignment_checker=None,
        gzipped: bool = False,
    ) -> None:
        super().__init__(output_file_name, params, assignment_checker)
        self.gzipped = gzipped
        self._closed = False
        if gzipped:
            self.output_file_name = output_file_name + ".gz"
            self.output_file = gzip.open(self.output_file_name, "wt")
        else:
            self.output_file = open(self.output_file_name, "w")

    def close(self) -> None:
        if not self._closed:
            self._closed = True
            f = getattr(self, "output_file", None)
            if f is not None and not f.closed:
                f.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def flush(self) -> None:
        if not self.gzipped and not self._closed and not self.output_file.closed:
            self.output_file.flush()


class ReadAssignmentCompositePrinter:
    def __init__(self, printers: Sequence) -> None:
        self.printers = list(printers)

    def __enter__(self) -> "ReadAssignmentCompositePrinter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def add_read_info(self, read_assignment) -> None:
        for p in self.printers:
            p.add_read_info(read_assignment)

    def flush(self) -> None:
        for p in self.printers:
            p.flush()

    def close(self) -> None:
        for p in self.printers:
            if hasattr(p, "close"):
                p.close()


class VoidPrinter:
    def __enter__(self) -> "VoidPrinter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def add_read_info(self, _) -> None:
        pass

    def flush(self) -> None:
        pass

    def close(self) -> None:
        pass


VoidAssignmentPrinter = VoidPrinter


class BEDPrinter(TextFileAssignmentPrinter):
    def __init__(
        self,
        output_file_name: str,
        params,
        print_corrected: bool = False,
        assignment_checker=None,
        gzipped: bool = False,
    ) -> None:
        super().__init__(output_file_name, params, assignment_checker, gzipped=gzipped)
        self.print_corrected = print_corrected
        self.output_file.write(
            "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\n"
        )

    def add_read_info(self, read_assignment) -> None:
        if (
            read_assignment is None
            or read_assignment.assignment_type is None
            or not hasattr(read_assignment, "gene_info")
            or read_assignment.gene_info is None
        ):
            return
        if self.assignment_checker is None or not self.assignment_checker.check(read_assignment):
            return

        strand = read_assignment.mapped_strand
        chr_id = read_assignment.gene_info.chr_id
        exon_blocks = read_assignment.corrected_exons if self.print_corrected else read_assignment.exons
        self.output_file.write(
            "%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n"
            % (
                chr_id,
                exon_blocks[0][0] - 1,
                exon_blocks[-1][1],
                read_assignment.read_id,
                strand,
                exon_blocks[0][0] - 1,
                exon_blocks[0][0] - 1,
                0,
                len(exon_blocks),
                ",".join(str(e[1] - e[0] + 1) for e in exon_blocks),
                ",".join(str(e[0] - exon_blocks[0][0]) for e in exon_blocks),
            )
        )


class TmpFileAssignmentPrinter(AbstractAssignmentPrinter):
    GENE_INFO = 255
    READ_ASSIGNMENT = 255 << 8

    def __init__(self, output_file_name: str, params) -> None:
        super().__init__(output_file_name, params)
        self.dumper = open(self.output_file_name, "wb")

    def close(self) -> None:
        dumper = getattr(self, "dumper", None)
        if dumper is not None and not dumper.closed:
            write_short_int(SHORT_TERMINATION_INT, dumper)
            dumper.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def add_gene_info(self, gene_info) -> None:
        write_short_int(self.GENE_INFO, self.dumper)
        gene_info.serialize(self.dumper)

    def add_read_info(self, read_assignment) -> None:
        if read_assignment.assignment_type is None or read_assignment.isoform_matches is None:
            return
        write_short_int(self.READ_ASSIGNMENT, self.dumper)
        read_assignment.serialize(self.dumper)

    def flush(self) -> None:
        if self.dumper is not None and not self.dumper.closed:
            self.dumper.flush()


class BaseTmpFileAssignmentLoader:
    def __init__(self, input_file_name: str) -> None:
        self.loader = open(input_file_name, "rb")
        self.current_id: Optional[int] = None
        self._read_id()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def close(self) -> None:
        loader = getattr(self, "loader", None)
        if loader is not None and not loader.closed:
            loader.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _read_id(self) -> None:
        self.current_id = read_short_int(self.loader)

    def has_next(self) -> bool:
        return self.current_id != SHORT_TERMINATION_INT

    def is_gene_info(self) -> bool:
        return self.current_id == TmpFileAssignmentPrinter.GENE_INFO

    def is_read_assignment(self) -> bool:
        return self.current_id == TmpFileAssignmentPrinter.READ_ASSIGNMENT


class NormalTmpFileAssignmentLoader(BaseTmpFileAssignmentLoader):
    def __init__(self, input_file_name: str, genedb, chr_record, string_pools=None) -> None:
        super().__init__(input_file_name)
        self.genedb = genedb
        self.chr_record = chr_record
        self.current_gene_info = None
        self.string_pools = string_pools

    def get_object(self):
        if self.is_gene_info():
            self.current_gene_info = GeneInfo.deserialize(self.loader, self.genedb)
            if self.chr_record is not None:
                self.current_gene_info.set_reference_sequence(
                    self.current_gene_info.all_read_region_start,
                    self.current_gene_info.all_read_region_end,
                    self.chr_record,
                )
            self._read_id()
            return self.current_gene_info
        if self.is_read_assignment():
            assert self.current_gene_info is not None
            assignment = ReadAssignment.deserialize(self.loader, self.current_gene_info, self.string_pools)
            self._read_id()
            return assignment
        return None


class GeneListTmpFileAssignmentLoader(BaseTmpFileAssignmentLoader):
    def __init__(self, input_file_name: str, string_pools=None) -> None:
        super().__init__(input_file_name)
        self.current_gene_info = None
        self.string_pools = string_pools

    def get_object(self):
        if self.is_gene_info():
            self.current_gene_info = GeneList.deserialize(self.loader)
            self._read_id()
            return self.current_gene_info
        if self.is_read_assignment():
            assert self.current_gene_info is not None
            assignment = ReadAssignment.deserialize(self.loader, self.current_gene_info, self.string_pools)
            self._read_id()
            return assignment
        return None


class QuickTmpFileAssignmentLoader(BaseTmpFileAssignmentLoader):
    def __init__(self, input_file_name: str, string_pools=None) -> None:
        super().__init__(input_file_name)
        self.string_pools = string_pools

    def get_object(self):
        if self.is_gene_info():
            GeneInfo.deserialize(self.loader, None)
            self._read_id()
            return None
        if self.is_read_assignment():
            assignment = BasicReadAssignment.deserialize_from_read_assignment(self.loader, self.string_pools)
            self._read_id()
            return assignment
        return None


class BasicTSVAssignmentPrinter(TextFileAssignmentPrinter):
    def __init__(
        self,
        output_file_name: str,
        params,
        io_support: "IOSupport",
        additional_header: str = "",
        gzipped: bool = False,
    ) -> None:
        super().__init__(output_file_name, params, gzipped=gzipped)
        self.header = (
            "read_id\tchr\tstrand\tisoform_id\tgene_id"
            "\tassignment_type\tassignment_events\texons\tadditional_info\tgroups\n"
        )
        self.output_file.write(additional_header)
        self.output_file.write(self.header)
        self.io_support = io_support

    @staticmethod
    def unmatched_line(read_assignment: ReadAssignment, additional_info: List[str]) -> str:
        read_exons = read_assignment.exons
        strand = read_assignment.strand
        line = (
            f"{read_assignment.read_id}\t{read_assignment.chr_id}\t{strand}\t.\t.\t"
            f"{read_assignment.assignment_type.name}\t.\t{range_list_to_str(read_exons)}"
        )
        for attr, val in read_assignment.additional_attributes.items():
            additional_info.append(f"{attr}={val};")
        line += "\t" + " ".join(additional_info) if additional_info else "\t*"
        line += "\t%s\n" % ",".join(read_assignment.read_group)
        return line

    def add_read_info(self, read_assignment: Optional[ReadAssignment]) -> None:
        if read_assignment is None:
            return
        if self.assignment_checker is None or not self.assignment_checker.check(read_assignment):
            return
        if read_assignment.assignment_type is None:
            logger.warning("Empty assignment read id %s", read_assignment.read_id)
            return
        if read_assignment.exons is None:
            logger.warning(
                "Empty combined profile, read id %s, assignment type %s",
                read_assignment.read_id,
                str(read_assignment.assignment_type),
            )
            return

        read_exons = read_assignment.exons
        read_introns = junctions_from_blocks(read_exons)
        if not read_assignment.isoform_matches:
            self.output_file.write(BasicTSVAssignmentPrinter.unmatched_line(read_assignment, []))
            return

        for m in read_assignment.isoform_matches:
            if m.assigned_transcript is None:
                self.output_file.write(
                    BasicTSVAssignmentPrinter.unmatched_line(
                        read_assignment, [f"Classification={m.match_classification.name};"]
                    )
                )
                continue
            isoform_introns = read_assignment.gene_info.all_isoforms_introns[m.assigned_transcript]
            event_string = ",".join(
                match_subtype_to_str_with_additional_info(
                    x, read_assignment.strand, read_introns, isoform_introns
                )
                for x in m.match_subclassifications
            )
            strand = read_assignment.strand
            line = (
                f"{read_assignment.read_id}\t{read_assignment.chr_id}\t{strand}\t"
                f"{m.assigned_transcript}\t{m.assigned_gene}\t"
                f"{read_assignment.assignment_type.name}\t{event_string}\t{range_list_to_str(read_exons)}"
            )
            additional_info = [
                f"gene_assignment={read_assignment.gene_assignment_type.name};",
                f"PolyA={read_assignment.polyA_found};",
            ]
            if getattr(self.params, "cage", None) is not None:
                additional_info.append(f"CAGE={read_assignment.cage_found};")
            if getattr(self.params, "check_canonical", False) and read_assignment.gene_info.reference_region:
                if not read_introns:
                    additional_info.append("Canonical=Unspliced;")
                else:
                    all_canonical = self.io_support.check_sites_are_canonical(
                        read_introns, read_assignment.gene_info, strand
                    )
                    additional_info.append(f"Canonical={all_canonical};")
            additional_info.append(f"Classification={m.match_classification.name};")
            for attr, val in read_assignment.additional_attributes.items():
                additional_info.append(f"{attr}={val};")
            line += "\t" + " ".join(additional_info) if additional_info else "\t*"
            line += "\t%s\n" % ",".join(read_assignment.read_group)
            self.output_file.write(line)


class ReadInfoPrinter(TextFileAssignmentPrinter):
    """Unified read_info.tsv formatter with promoted columns and corrected exon coordinates."""

    HEADER = (
        "read_id\tchr\tstrand\tgene_id\tgene_assignment_type\tisoform_id\tisoform_assignment_type"
        "\tassignment_events\tclassification\texons\tpolyA\tCAGE\tcanonical"
        "\tbarcode\tumi\tcell_type\tgroups\tadditional\n"
    )

    def __init__(
        self,
        output_file_name: str,
        params,
        io_support: "IOSupport",
        additional_header: str = "",
        gzipped: bool = False,
    ) -> None:
        super().__init__(output_file_name, params, gzipped=gzipped)
        self.output_file.write(additional_header)
        self.output_file.write(self.HEADER)
        self.io_support = io_support

    def _get_exons_str(self, read_assignment: ReadAssignment) -> str:
        exon_blocks = read_assignment.corrected_exons if read_assignment.corrected_exons else read_assignment.exons
        return "." if exon_blocks is None else range_list_to_str(exon_blocks)

    def _get_canonical(self, read_introns, read_assignment: ReadAssignment, strand: str) -> str:
        if (
            not getattr(self.params, "check_canonical", False)
            or not read_assignment.gene_info
            or not read_assignment.gene_info.reference_region
        ):
            return "."
        if not read_introns:
            return "Unspliced"
        return str(self.io_support.check_sites_are_canonical(read_introns, read_assignment.gene_info, strand))

    def _get_additional(self, read_assignment: ReadAssignment, exclude_keys: set[str]) -> str:
        parts = [f"{k}={v};" for k, v in read_assignment.additional_attributes.items() if k not in exclude_keys]
        return " ".join(parts) if parts else "*"

    def _write_unmatched(
        self, read_assignment: ReadAssignment, classification: str = ".", gene_id: str = "."
    ) -> None:
        read_introns = junctions_from_blocks(read_assignment.exons) if read_assignment.exons else []
        strand = read_assignment.strand
        fields = [
            read_assignment.read_id,
            read_assignment.chr_id,
            strand,
            gene_id,
            read_assignment.gene_assignment_type.name,
            ".",
            read_assignment.assignment_type.name,
            ".",
            classification,
            self._get_exons_str(read_assignment),
            str(read_assignment.polyA_found),
            str(read_assignment.cage_found) if getattr(self.params, "cage", None) is not None else ".",
            self._get_canonical(read_introns, read_assignment, strand),
            read_assignment.barcode or ".",
            read_assignment.umi or ".",
            read_assignment.additional_attributes.get("cell_type", "."),
            ",".join(read_assignment.read_group) if read_assignment.read_group else ".",
            self._get_additional(read_assignment, {"cell_type"}),
        ]
        self.output_file.write("\t".join(fields) + "\n")

    def add_read_info(self, read_assignment: Optional[ReadAssignment]) -> None:
        if read_assignment is None:
            return
        if self.assignment_checker is None or not self.assignment_checker.check(read_assignment):
            return
        if read_assignment.assignment_type is None:
            logger.warning("Empty assignment read id %s", read_assignment.read_id)
            return
        if read_assignment.exons is None:
            logger.warning(
                "Empty combined profile, read id %s, assignment type %s",
                read_assignment.read_id,
                str(read_assignment.assignment_type),
            )
            return

        read_exons = read_assignment.exons
        read_introns = junctions_from_blocks(read_exons)
        if not read_assignment.isoform_matches:
            self._write_unmatched(read_assignment)
            return

        for m in read_assignment.isoform_matches:
            if m.assigned_transcript is None:
                self._write_unmatched(
                    read_assignment,
                    m.match_classification.name,
                    m.assigned_gene if m.assigned_gene else ".",
                )
                continue

            isoform_introns = read_assignment.gene_info.all_isoforms_introns[m.assigned_transcript]
            event_string = ",".join(
                match_subtype_to_str_with_additional_info(
                    x, read_assignment.strand, read_introns, isoform_introns
                )
                for x in m.match_subclassifications
            )
            strand = read_assignment.strand
            fields = [
                read_assignment.read_id,
                read_assignment.chr_id,
                strand,
                m.assigned_gene,
                read_assignment.gene_assignment_type.name,
                m.assigned_transcript,
                read_assignment.assignment_type.name,
                event_string,
                m.match_classification.name,
                self._get_exons_str(read_assignment),
                str(read_assignment.polyA_found),
                str(read_assignment.cage_found) if getattr(self.params, "cage", None) is not None else ".",
                self._get_canonical(read_introns, read_assignment, strand),
                read_assignment.barcode or ".",
                read_assignment.umi or ".",
                read_assignment.additional_attributes.get("cell_type", "."),
                ",".join(read_assignment.read_group) if read_assignment.read_group else ".",
                self._get_additional(read_assignment, {"cell_type"}),
            ]
            self.output_file.write("\t".join(fields) + "\n")


class SqantiTSVPrinter(TextFileAssignmentPrinter):
    def __init__(self, output_file_name: str, params, io_support: "IOSupport") -> None:
        super().__init__(output_file_name, params)
        self.header = (
            "isoform\tchrom\tstrand\tlength\texons\tstructural_category"
            "\tassociated_gene\tassociated_transcript\tref_length\tref_exons\tdiff_to_TSS\tdiff_to_TTS"
            "\tdiff_to_gene_TSS\tdiff_to_gene_TTS\tsubcategory\tRTS_stage\tall_canonical"
            "\tmin_sample_cov\tmin_cov\tmin_cov_pos\tsd_cov\tFL"
            "\tn_indels\tn_indels_junc\tbite\tiso_exp\tgene_exp\tratio_exp\tFSM_class"
            "\tcoding\tORF_length\tCDS_length\tCDS_start\tCDS_end"
            "\tCDS_genomic_start\tCDS_genomic_end\tpredicted_NMD"
            "\tperc_A_downstream_TTS\tseq_A_downstream_TTS\tdist_to_CAGE_peak"
            "\twithin_CAGE_peak\tdist_to_polyA_site\twithin_polyA_site\tpolyA_motif\tpolyA_dist"
            "\tpolyA_motif_found\tORF_seq\tratio_TSS\n"
        )
        self.output_file.write(self.header)
        self.output_file.flush()
        self.io_support = io_support

    def add_read_info(self, read_assignment: Optional[ReadAssignment]) -> None:
        if read_assignment is None or read_assignment.assignment_type == ReadAssignmentType.noninformative:
            return

        gene_info = read_assignment.gene_info
        if not read_assignment.isoform_matches:
            assigned = False
            gene_id = "NA"
            transcript_id = "NA"
            sqanti_classification = MatchClassification.intergenic
            match = None
        else:
            match = read_assignment.isoform_matches[0]
            gene_id = match.assigned_gene
            transcript_id = match.assigned_transcript
            sqanti_classification = match.match_classification
            if transcript_id is None:
                assigned = False
                transcript_id = "NA"
                gene_id = "NA"
                sqanti_classification = MatchClassification.intergenic
            else:
                assigned = True

        strand = read_assignment.strand
        if assigned:
            dist_to_tss = self.io_support.count_tss_dist(read_assignment, transcript_id)
            dist_to_tts = self.io_support.count_tts_dist(read_assignment, transcript_id)
            dist_to_gene_tss, dist_to_gene_tts = self.io_support.find_closests_tsts(read_assignment)
        else:
            dist_to_tss = dist_to_tts = "NA"
            dist_to_gene_tss = dist_to_gene_tts = "NA"
        if strand == "-" and assigned:
            dist_to_tss, dist_to_tts = dist_to_tts, dist_to_tss
            dist_to_gene_tss, dist_to_gene_tts = dist_to_gene_tts, dist_to_gene_tss

        subtypes = (
            ";".join(x.event_type.name for x in match.match_subclassifications) or "."
            if assigned and match is not None
            else "."
        )
        indel_count = read_assignment.additional_info.get("indel_count", "-1")
        junctions_with_indels = read_assignment.additional_info.get("junctions_with_indels", "-1")
        bite = self.io_support.check_all_sites_match_reference(read_assignment)
        ref_cds_start, ref_cds_end = (
            self.io_support.find_ref_CDS_region(gene_info, transcript_id) if assigned else ("NA", "NA")
        )

        if gene_info.reference_region:
            read_introns = junctions_from_blocks(read_assignment.exons)
            all_canonical = str(self.io_support.check_sites_are_canonical(read_introns, gene_info, strand))
            seq_a_down, perc_a_down = self.io_support.check_downstream_polya(
                (read_assignment.start(), read_assignment.end()), gene_info, strand
            )
            perc_a_str = f"{perc_a_down:0.2f}"
        else:
            all_canonical = "NA"
            perc_a_str = "NA"
            seq_a_down = "NA"

        transcript_len = "NA" if not assigned else gene_info.total_transcript_length(transcript_id)
        transcript_exon_count = "NA" if not assigned else gene_info.transcript_exon_count(transcript_id)
        fsm_class = read_assignment.additional_info["FSM_class"]

        values = [
            read_assignment.read_id, gene_info.chr_id, strand, read_assignment.length(),
            read_assignment.exon_count(), sqanti_classification.name, gene_id, transcript_id,
            transcript_len, transcript_exon_count,
            dist_to_tss, dist_to_tts, dist_to_gene_tss, dist_to_gene_tts, subtypes, "FALSE",
            all_canonical, "NA", "NA", "NA", "NA", "NA",
            indel_count, junctions_with_indels, bite, "NA", "NA", "NA", fsm_class,
            "NA", "NA", "NA", "NA", "NA",
            ref_cds_start, ref_cds_end, "NA", perc_a_str, seq_a_down,
            "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
        ]
        self.output_file.write("\t".join(str(x) for x in values) + "\n")


class IOSupport:
    """SQANTI / transcript-distance / splice-site / downstream polyA helper supporting both .iqidx and gffutils."""

    def __init__(self, params) -> None:
        self.params = params

    def count_tss_dist(self, read_assignment: ReadAssignment, transcript_id: str) -> int:
        transcript_start = read_assignment.gene_info.transcript_start(transcript_id)
        read_start = read_assignment.start()
        if read_start < transcript_start:
            return -sum_intervals_to_point(read_assignment.exons, transcript_start)
        return sum_intervals_to_point(read_assignment.gene_info.all_isoforms_exons[transcript_id], read_start)

    def count_tts_dist(self, read_assignment: ReadAssignment, transcript_id: str) -> int:
        transcript_end = read_assignment.gene_info.transcript_end(transcript_id)
        read_end = read_assignment.end()
        if read_end > transcript_end:
            return -sum_intervals_from_point(read_assignment.exons, transcript_end)
        return sum_intervals_from_point(read_assignment.gene_info.all_isoforms_exons[transcript_id], read_end)

    def find_closests_tsts(self, read_assignment: ReadAssignment) -> Tuple[int, int]:
        gene_info = read_assignment.gene_info
        gene_id = read_assignment.isoform_matches[0].assigned_gene

        # Prefer in-memory GeneInfo map populated directly from .iqidx without SQLite queries.
        tx_ids = [tid for tid, gid in gene_info.gene_id_map.items() if gid == gene_id]
        if tx_ids:
            tx_ids.sort(key=lambda tid: (gene_info.transcript_start(tid), tid))
        elif gene_info.db is not None and hasattr(gene_info.db, "children"):
            tx_ids = [
                t.id
                for t in gene_info.db.children(gene_info.db[gene_id], featuretype=("transcript", "mRNA"), order_by="start")
            ]

        dists_tss = [self.count_tss_dist(read_assignment, tid) for tid in tx_ids]
        dists_tts = [self.count_tts_dist(read_assignment, tid) for tid in tx_ids]
        return (
            dists_tss[argmin([abs(x) for x in dists_tss])],
            dists_tts[argmin([abs(x) for x in dists_tts])],
        )

    def find_ref_CDS_region(self, gene_info, transcript_id: str) -> Tuple[int, int]:
        other = gene_info.other_features.get(transcript_id)
        cds_intervals = sorted((s, e) for s, e, ftype in other if ftype == "CDS") if other else []
        if not cds_intervals and gene_info.db is not None and hasattr(gene_info.db, "children"):
            cds_intervals = [
                (c.start, c.end)
                for c in gene_info.db.children(gene_info.db[transcript_id], featuretype="CDS", order_by="start")
            ]
        if not cds_intervals:
            return -1, -1
        if gene_info.isoform_strands[transcript_id] == "-":
            return cds_intervals[-1][1], cds_intervals[0][0]
        return cds_intervals[0][0], cds_intervals[-1][1]

    def check_all_sites_match_reference(self, read_assignment: ReadAssignment) -> str:
        if len(read_assignment.exons) == 1:
            return "Unspliced"
        return str(read_assignment.introns_match)

    def check_sites_are_canonical(
        self, read_introns: Iterable[Tuple[int, int]], gene_info, strand: str
    ) -> bool:
        for intron in read_introns:
            if intron not in gene_info.canonical_sites:
                intron_left_pos = intron[0] - gene_info.all_read_region_start
                intron_right_pos = intron[1] - gene_info.all_read_region_start
                left_site = gene_info.reference_region[intron_left_pos : intron_left_pos + 2]
                right_site = gene_info.reference_region[intron_right_pos - 1 : intron_right_pos + 1]
                if strand == "+":
                    gene_info.canonical_sites[intron] = (left_site, right_site) in CANONICAL_FWD_SITES
                else:
                    gene_info.canonical_sites[intron] = (left_site, right_site) in CANONICAL_REV_SITES
            if not gene_info.canonical_sites[intron]:
                return False
        return True

    def add_canonical_info_for_model(self, model, gene_info) -> None:
        if not gene_info.reference_region:
            return
        key_word = "Canonical"
        if model.check_additional(key_word):
            return
        model_introns = junctions_from_blocks(model.exon_blocks)
        if not model_introns:
            canonical_info = "Unspliced"
        else:
            canonical_info = str(self.check_sites_are_canonical(model_introns, gene_info, model.strand))
        model.add_additional_attribute(key_word, canonical_info)

    def add_canonical_info(self, model_storage, gene_info) -> None:
        for m in model_storage:
            self.add_canonical_info_for_model(m, gene_info)

    def check_downstream_polya(
        self, read_coords: Tuple[int, int], gene_info, strand: str
    ) -> Tuple[str, float]:
        up_len = self.params.upstream_region_len
        if strand == "+":
            read_end = read_coords[1] - gene_info.all_read_region_start + 1
            seq = gene_info.reference_region[read_end : read_end + up_len]
            a_pct = float(seq.upper().count("A")) / float(up_len)
        else:
            read_start = read_coords[0] - gene_info.all_read_region_start
            seq = gene_info.reference_region[read_start - up_len : read_start]
            a_pct = float(seq.upper().count("T")) / float(up_len)
        return seq, a_pct
