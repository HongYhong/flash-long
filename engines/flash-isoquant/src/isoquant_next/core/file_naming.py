"""Canonical file naming conventions compatible with IsoQuant 4.0.0."""

from __future__ import annotations

import os
from typing import Callable, Iterable, Optional


def convert_chr_id_to_file_name_str(chr_id: str) -> str:
    return chr_id.replace("/", "_")


def reads_collected_lock_file_name(sample_out_raw: str, chr_id: str) -> str:
    return f"{sample_out_raw}_{convert_chr_id_to_file_name_str(chr_id)}_collected"


def reads_processed_lock_file_name(dump_filename: str, chr_id: str) -> str:
    return f"{dump_filename}_{convert_chr_id_to_file_name_str(chr_id)}_processed"


def read_group_lock_filename(sample) -> str:
    return f"{sample.read_group_file}_lock"


def split_barcodes_lock_filename(sample) -> str:
    return f"{sample.barcodes_split_reads}_lock"


def clean_locks(chr_ids: Iterable[str], base_name: str, fname_function: Callable[[str, str], str]) -> None:
    for chr_id in chr_ids:
        fname = fname_function(base_name, chr_id)
        if os.path.exists(fname):
            os.remove(fname)


def saves_file_name(out_raw_file: str, chr_id: str) -> str:
    return f"{out_raw_file}_{convert_chr_id_to_file_name_str(chr_id)}"


def multimappers_file_name(out_raw_file: str, chr_id: str) -> str:
    return f"{out_raw_file}_multimappers_{convert_chr_id_to_file_name_str(chr_id)}"


def filtered_reads_file_name(out_raw_file: str, chr_id: str) -> str:
    return f"{out_raw_file}_filtered_{convert_chr_id_to_file_name_str(chr_id)}"


def umi_filtered_reads_file_name(out_umi_filtered_tmp: str, chr_id: str, edit_distance: int) -> str:
    return f"{out_umi_filtered_tmp}_{convert_chr_id_to_file_name_str(chr_id)}_ED{edit_distance}"


def umi_filtered_lock_file_name(out_umi_filtered_done: str, chr_id: str, edit_distance: int) -> str:
    return f"{out_umi_filtered_done}_{convert_chr_id_to_file_name_str(chr_id)}_ED{edit_distance}"


def allinfo_file_name(out_umi_filtered_tmp: str, chr_id: str, edit_distance: int) -> str:
    return umi_filtered_reads_file_name(out_umi_filtered_tmp, chr_id, edit_distance) + ".allinfo"


def allinfo_stats_file_name(out_umi_filtered_tmp: str, chr_id: str, edit_distance: int) -> str:
    return umi_filtered_reads_file_name(out_umi_filtered_tmp, chr_id, edit_distance) + ".stats.tsv"


def umi_filtered_global_lock_file_name(out_umi_filtered_done: str) -> str:
    return f"{out_umi_filtered_done}.lock"


def umi_barcode2barcode_prefix(out_umi_filtered_done: str, col_index: int) -> str:
    return f"{out_umi_filtered_done}.barcode_barcode_col{col_index}"


def umi_barcode2barcode_global_lock(out_umi_filtered_done: str, col_index: int) -> str:
    return umi_barcode2barcode_prefix(out_umi_filtered_done, col_index) + ".lock"


def dynamic_pools_file_name(out_raw_file: str, chr_id: str) -> str:
    return f"{out_raw_file}_dynamic_pools_{convert_chr_id_to_file_name_str(chr_id)}"


def barcode_pools_file_name(out_raw_file: str, chr_id: str) -> str:
    return f"{out_raw_file}_barcode_pools_{convert_chr_id_to_file_name_str(chr_id)}"


def read_groups_file_name(save_file: str) -> str:
    return f"{save_file}_groups"


def bamstat_file_name(save_file: str) -> str:
    return f"{save_file}_bamstat"


def info_file_name(out_raw_file: str) -> str:
    return f"{out_raw_file}_info"


def collection_lock_file_name(out_raw_file: str) -> str:
    return f"{out_raw_file}_lock"


def read_stat_file_name(saves_file: str) -> str:
    return f"{saves_file}_read_stat"


def transcript_stat_file_name(saves_file: str) -> str:
    return f"{saves_file}_transcript_stat"


def read_group_file_from_saves(saves_prefix: str) -> Optional[str]:
    if saves_prefix.endswith(".save"):
        return saves_prefix[:-5] + ".read_group"
    return None


def umi_output_prefix(out_umi_filtered: str, edit_distance: int) -> str:
    if edit_distance < 0:
        return f"{out_umi_filtered}.ALL"
    return f"{out_umi_filtered}.ED{edit_distance}"


def counts_prefix(output_prefix: str) -> str:
    return f"{output_prefix}_counts"


def tpm_prefix(output_prefix: str) -> str:
    return f"{output_prefix}_tpm"


def counts_file_name(counts_prefix_str: str, linear: bool = False) -> str:
    return counts_prefix_str + (".linear.tsv" if linear else ".tsv")


def tpm_file_name(tpm_prefix_str: str) -> str:
    return f"{tpm_prefix_str}.tsv"


def counts_stats_file_name(counts_file: str) -> str:
    return f"{counts_file}.stats"


def counts_usable_file_name(counts_file: str) -> str:
    return f"{counts_file}.usable"


def mtx_matrix_file(output_prefix: str) -> str:
    return f"{output_prefix}.matrix.mtx"


def mtx_features_file(output_prefix: str) -> str:
    return f"{output_prefix}.features.tsv"


def mtx_barcodes_file(output_prefix: str) -> str:
    return f"{output_prefix}.barcodes.tsv"


def mtx_include_matrix_file(output_prefix: str) -> str:
    return f"{output_prefix}.include.matrix.mtx"


def mtx_exclude_matrix_file(output_prefix: str) -> str:
    return f"{output_prefix}.exclude.matrix.mtx"
