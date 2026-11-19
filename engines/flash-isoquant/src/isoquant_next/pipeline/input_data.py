"""Sample and input dataset containers compatible with FASTQ, BAM, unmapped BAM, YAML, and saved runs."""

from __future__ import annotations

from collections import defaultdict
from enum import Enum, unique
import logging
import os
from typing import Dict, List, Optional

import yaml

from isoquant_next.core.error_codes import IsoQuantError, IsoQuantExitCode
from isoquant_next.core.file_naming import convert_chr_id_to_file_name_str
from isoquant_next.core.file_utils import normalize_path

logger = logging.getLogger("IsoQuant")


@unique
class InputDataType(Enum):
    undefined = 0
    fastq = 1
    bam = 2
    unmapped_bam = 3
    save = 10

    def needs_mapping(self) -> bool:
        return self in (InputDataType.fastq, InputDataType.unmapped_bam)


def check_input_type(fname: str, input_type: InputDataType) -> None:
    basename_plus_inner_ext, outer_ext = os.path.splitext(fname.lower())
    if outer_ext not in (".zip", ".gz", ".gzip", ".bz2", ".bzip2"):
        basename_plus_inner_ext = fname
    _base, inner_ext = os.path.splitext(basename_plus_inner_ext)
    inner_ext = inner_ext.lower()
    if inner_ext in (".fastq", ".fasta", ".fa", ".fq", ".fna"):
        if input_type != InputDataType.fastq:
            raise IsoQuantError(
                IsoQuantExitCode.INVALID_FILE_FORMAT,
                f"Wrong file extension detected {fname}; use FASTQ/FASTA only with --fastq.",
            )
    elif inner_ext == ".bam":
        if input_type not in (InputDataType.bam, InputDataType.unmapped_bam):
            raise IsoQuantError(
                IsoQuantExitCode.INVALID_FILE_FORMAT,
                f"Wrong file extension detected {fname}; use BAM only with --bam/--unmapped_bam.",
            )
    else:
        raise IsoQuantError(
            IsoQuantExitCode.INVALID_FILE_FORMAT,
            f"Unsupported input format '{inner_ext}' for {fname}; supported: FASTQ, FASTA, BAM.",
        )


class SampleData:
    """Pure-data per-sample input specification and deterministic output/aux paths."""

    def __init__(
        self,
        file_list: List[List[str]],
        prefix: str,
        out_dir: str,
        readable_names_dict: Dict[str, str],
        illumina_bam: Optional[List[str]] = None,
        barcoded_reads: Optional[List[str]] = None,
    ) -> None:
        self.file_list = file_list
        self.readable_names_dict = dict(readable_names_dict)
        self.illumina_bam = illumina_bam
        self.prefix = prefix
        self.out_dir = out_dir
        self.aux_dir = os.path.join(self.out_dir, "aux")
        self.barcoded_reads = list(barcoded_reads) if barcoded_reads else []
        self.use_technical_replicas = False
        self._init_paths()

    def _make_path(self, name: str) -> str:
        return os.path.join(self.out_dir, name)

    def _make_aux_path(self, name: str) -> str:
        return os.path.join(self.aux_dir, name)

    def _init_paths(self) -> None:
        p = self.prefix
        self.out_assigned_tsv = self._make_path(f"{p}.read_assignments.tsv")
        self.out_assigned_tsv_result = self.out_assigned_tsv
        self.out_read_info_tsv = self._make_path(f"{p}.read_info.tsv")
        self.out_transcript_model_reads_tsv = self._make_path(f"{p}.transcript_model_reads.tsv")
        self.out_raw_file = self._make_aux_path(f"{p}.save")
        self.read_group_file = self._make_aux_path(f"{p}.read_group")
        self.out_corrected_bed = self._make_path(f"{p}.corrected_reads.bed")
        self.out_alt_tsv = self._make_path(f"{p}.read_assignments.SQANTI-like.tsv")
        self.out_gene_counts_tsv = self._make_path(f"{p}.gene")
        self.out_transcript_counts_tsv = self._make_path(f"{p}.transcript")
        self.out_transcript_model_counts_tsv = self._make_path(f"{p}.discovered_transcript")
        self.out_gene_model_counts_tsv = self._make_path(f"{p}.discovered_gene")
        self.out_transcript_model_grouped_counts_tsv = self._make_path(f"{p}.discovered_transcript_grouped")
        self.out_gene_model_grouped_counts_tsv = self._make_path(f"{p}.discovered_gene_grouped")
        self.out_exon_counts_tsv = self._make_path(f"{p}.exon")
        self.out_exon_splice_site_counts_tsv = self._make_path(f"{p}.exon_splice_site")
        self.out_old_exon_counts_tsv = self._make_path(f"{p}.old_exon")
        self.out_intron_counts_tsv = self._make_path(f"{p}.splice_junction")
        self.out_intron_retention_counts_tsv = self._make_path(f"{p}.intron_retention")
        self.out_intron_retention_grouped_counts_tsv = self._make_path(f"{p}.intron_retention_grouped")
        self.out_gene_grouped_counts_tsv = self._make_path(f"{p}.gene_grouped")
        self.out_transcript_grouped_counts_tsv = self._make_path(f"{p}.transcript_grouped")
        self.out_exon_grouped_counts_tsv = self._make_path(f"{p}.exon_grouped")
        self.out_exon_splice_site_grouped_counts_tsv = self._make_path(f"{p}.exon_splice_site_grouped")
        self.out_old_exon_grouped_counts_tsv = self._make_path(f"{p}.old_exon_grouped")
        self.out_intron_grouped_counts_tsv = self._make_path(f"{p}.splice_junction_grouped")
        self.out_t2t_tsv = self._make_path(f"{p}.novel_vs_known.SQANTI-like.tsv")
        self.out_polya_prediction_tsv = self._make_path(f"{p}.polyA_prediction.tsv")
        self.out_tss_prediction_tsv = self._make_path(f"{p}.TSS_prediction.tsv")
        self.out_polya_prediction_grouped_tsv = self._make_path(f"{p}.polyA_prediction_grouped")
        self.out_tss_prediction_grouped_tsv = self._make_path(f"{p}.TSS_prediction_grouped")
        self.barcodes_tsv = self._make_path(f"{p}.barcoded_reads")
        self.barcodes_done = self._make_aux_path(f"{p}.barcodes_done")
        self.barcodes_split_reads = self._make_aux_path(f"{p}.split_barcodes")
        self.out_umi_filtered = self._make_path(f"{p}.UMI_filtered")
        self.out_umi_filtered_tmp = self._make_aux_path(f"{p}.UMI_filtered")
        self.out_umi_filtered_done = self._make_aux_path(f"{p}.UMI_filtered.done")
        self.split_reads_fasta = self._make_path(f"{p}.split_reads")

    def get_save_file(self, chr_id: str) -> str:
        from isoquant_next.core.file_naming import saves_file_name

        return saves_file_name(self.out_raw_file, chr_id)

    def get_bamstat_file(self, chr_id: str) -> str:
        from isoquant_next.core.file_naming import bamstat_file_name

        return bamstat_file_name(self.get_save_file(chr_id))

    def get_filtered_reads_file(self, chr_id: str) -> str:
        from isoquant_next.core.file_naming import filtered_reads_file_name

        return filtered_reads_file_name(self.out_raw_file, chr_id)

    def get_dynamic_pools_file(self, chr_id: str) -> str:
        from isoquant_next.core.file_naming import dynamic_pools_file_name

        return dynamic_pools_file_name(self.out_raw_file, chr_id)

    def get_barcode_pools_file(self, chr_id: str) -> str:
        from isoquant_next.core.file_naming import barcode_pools_file_name

        return barcode_pools_file_name(self.out_raw_file, chr_id)

    def get_barcodes_split_file(self, chr_id: str) -> str:
        return f"{self.barcodes_split_reads}_{convert_chr_id_to_file_name_str(chr_id)}"

    def get_chr_prefix(self, chr_id: str) -> str:
        return f"{self.prefix}_{convert_chr_id_to_file_name_str(chr_id)}"

    def get_corrected_bed_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.corrected_reads.bed")

    def get_assigned_tsv_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.read_assignments.tsv")

    def get_read_info_tsv_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.read_info.tsv")

    def get_transcript_model_reads_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.transcript_model_reads.tsv")

    def get_gene_counts_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.gene")

    def get_transcript_counts_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.transcript")

    def get_exon_counts_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.exon")

    def get_old_exon_counts_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.old_exon")

    def get_exon_splice_site_counts_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.exon_splice_site")

    def get_intron_counts_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.splice_junction")

    def get_intron_retention_counts_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.intron_retention")

    def get_transcript_model_counts_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.discovered_transcript")

    def get_gene_model_counts_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.discovered_gene")

    def get_t2t_tsv_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.novel_vs_known.SQANTI-like.tsv")

    def get_polya_prediction_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.polyA_prediction.tsv")

    def get_tss_prediction_file(self, chr_id: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.TSS_prediction.tsv")

    def get_grouped_counts_file(self, chr_id: str, feature: str, strategy_name: str) -> str:
        return self._make_path(f"{self.get_chr_prefix(chr_id)}.{feature}_grouped_{strategy_name}")

    def get_info_file(self) -> str:
        from isoquant_next.core.file_naming import info_file_name

        return info_file_name(self.out_raw_file)

    def get_collection_lock_file(self) -> str:
        from isoquant_next.core.file_naming import collection_lock_file_name

        return collection_lock_file_name(self.out_raw_file)


class InputDataStorage:
    """Validates and normalizes CLI/YAML inputs into per-sample SampleData records."""

    def __init__(self, args) -> None:
        self.samples: List[SampleData] = []
        self.input_type = InputDataType.undefined
        self.experiment_prefix = args.prefix
        readable_names_dict: Dict[str, Dict[str, str]] = defaultdict(dict)
        sample_files: List[List[List[str]]] = []
        experiment_names: List[str] = []
        illumina_bam: List[Optional[List[str]]] = []

        if getattr(args, "fastq", None) is not None:
            self.input_type = InputDataType.fastq
            self._populate_single_experiment(
                args.fastq, args.labels, args.illumina_bam, sample_files, experiment_names, readable_names_dict, illumina_bam
            )
        elif getattr(args, "bam", None) is not None:
            self.input_type = InputDataType.bam
            self._populate_single_experiment(
                args.bam, args.labels, args.illumina_bam, sample_files, experiment_names, readable_names_dict, illumina_bam
            )
        elif getattr(args, "unmapped_bam", None) is not None:
            self.input_type = InputDataType.unmapped_bam
            self._populate_single_experiment(
                args.unmapped_bam, args.labels, args.illumina_bam, sample_files, experiment_names, readable_names_dict, illumina_bam
            )
        elif getattr(args, "read_assignments", None) is not None:
            self.input_type = InputDataType.save
            for i, save_file in enumerate(args.read_assignments):
                sample_files.append([[save_file]])
                experiment_names.append(f"{self.experiment_prefix}{i}")
                illumina_bam.append(None)
        elif getattr(args, "yaml", None) is not None:
            sample_files, experiment_names, readable_names_dict, illumina_bam = self._parse_yaml(args.yaml)
            if getattr(args, "labels", None):
                logger.warning("--labels option has no effect when files are provided via yaml file")
        else:
            raise IsoQuantError(IsoQuantExitCode.NO_INPUT_DATA, "Input data was not specified")

        for i, files in enumerate(sample_files):
            exp_name = experiment_names[i]
            self.samples.append(
                SampleData(
                    files,
                    exp_name,
                    os.path.join(args.output, exp_name),
                    readable_names_dict[exp_name],
                    illumina_bam[i],
                )
            )

    def _populate_single_experiment(
        self,
        paths: List[str],
        labels: Optional[List[str]],
        illumina: Optional[List[str]],
        sample_files: List[List[List[str]]],
        experiment_names: List[str],
        readable_names_dict: Dict[str, Dict[str, str]],
        illumina_bam: List[Optional[List[str]]],
    ) -> None:
        exp_name = self.experiment_prefix
        if labels and len(labels) != len(paths):
            raise IsoQuantError(
                IsoQuantExitCode.LABEL_COUNT_MISMATCH,
                "Number of labels is not equal to the number of input files",
            )
        entries: List[List[str]] = []
        for i, p in enumerate(paths):
            check_input_type(p, self.input_type)
            if p in readable_names_dict[exp_name]:
                raise IsoQuantError(
                    IsoQuantExitCode.DUPLICATE_FILES,
                    f"File {p} is used multiple times in a single experiment",
                )
            entries.append([p])
            readable_names_dict[exp_name][p] = (
                labels[i] if labels else os.path.splitext(os.path.basename(p))[0]
            )
        sample_files.append(entries)
        experiment_names.append(exp_name)
        illumina_bam.append(illumina)

    def _parse_yaml(self, yaml_file_path: str):
        with open(yaml_file_path, "r", encoding="utf-8") as fh:
            con = yaml.safe_load(fh)
        if not con or not isinstance(con, list) or "data format" not in con[0]:
            raise IsoQuantError(
                IsoQuantExitCode.YAML_PARSING_ERROR,
                "YAML first entry must specify 'data format' (fastq, fasta, bam, or unmapped_bam)",
            )
        fmt = con[0]["data format"]
        fmt_map = {
            "bam": InputDataType.bam,
            "unmapped_bam": InputDataType.unmapped_bam,
            "fastq": InputDataType.fastq,
            "fasta": InputDataType.fastq,
        }
        if fmt not in fmt_map:
            raise IsoQuantError(
                IsoQuantExitCode.INVALID_FILE_FORMAT,
                f"Unsupported YAML data format: {fmt}",
            )
        self.input_type = fmt_map[fmt]

        sample_files: List[List[List[str]]] = []
        experiment_names: List[str] = []
        illumina_bam: List[Optional[List[str]]] = []
        readable_names_dict: Dict[str, Dict[str, str]] = defaultdict(dict)

        for idx, sample in enumerate(con[1:]):
            name = sample.get("name") or f"{self.experiment_prefix}{idx}"
            if name in experiment_names:
                alt = f"{self.experiment_prefix}{idx}"
                if alt == name:
                    raise IsoQuantError(IsoQuantExitCode.DUPLICATE_FILES, f"Duplicate experiment name: {name}")
                name = alt
            if "long read files" not in sample:
                raise IsoQuantError(IsoQuantExitCode.YAML_PARSING_ERROR, f"Experiment {name} missing 'long read files'")
            files = [normalize_path(yaml_file_path, f) for f in sample["long read files"]]
            labels = sample.get("labels")
            if labels is not None and len(labels) != len(files):
                raise IsoQuantError(IsoQuantExitCode.LABEL_COUNT_MISMATCH, "YAML labels length != long read files length")
            for fi, fname in enumerate(files):
                if fname in readable_names_dict[name]:
                    raise IsoQuantError(IsoQuantExitCode.DUPLICATE_FILES, f"Duplicate file {fname} in {name}")
                check_input_type(fname, self.input_type)
                readable_names_dict[name][fname] = (
                    labels[fi] if labels is not None else os.path.splitext(os.path.basename(fname))[0]
                )
            sample_files.append([[f] for f in files])
            experiment_names.append(name)
            ib = sample.get("illumina bam")
            illumina_bam.append([normalize_path(yaml_file_path, x) for x in ib] if ib else None)

        return sample_files, experiment_names, readable_names_dict, illumina_bam

    def has_replicas(self) -> bool:
        return any(len(s.file_list) > 1 for s in self.samples)
