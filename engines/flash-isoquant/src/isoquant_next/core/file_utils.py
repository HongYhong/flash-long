"""File merging, path normalization, and resource-safe count merging helpers."""

from __future__ import annotations

import os
import re
import shutil
from collections import defaultdict
from typing import Iterable, List, Optional, TextIO

from isoquant_next.core.file_naming import convert_chr_id_to_file_name_str
from isoquant_next.core.utils import rreplace


def _natural_sort_key(path_str: str) -> list[int | str]:
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", path_str)]


def merge_file_list(fname: str, label: str, chr_ids: Iterable[str]) -> List[str]:
    """Build per-chromosome shard filenames without matching short labels inside basename suffixes."""
    directory, base = os.path.split(fname)
    if base.startswith(label):
        rest = base[len(label):]
        return [
            os.path.join(directory, f"{label}_{convert_chr_id_to_file_name_str(chr_id)}{rest}")
            for chr_id in chr_ids
        ]
    return [rreplace(fname, label, f"{label}_{convert_chr_id_to_file_name_str(chr_id)}") for chr_id in chr_ids]


def merge_files(
    file_name: str,
    label: str,
    chr_ids: Iterable[str],
    merged_file_handler: TextIO,
    copy_header: bool = True,
    header_lines: Optional[int] = None,
) -> None:
    file_names = sorted(merge_file_list(file_name, label, chr_ids), key=_natural_sort_key)
    for i, chr_file in enumerate(file_names):
        if not os.path.exists(chr_file):
            continue
        if header_lines is not None:
            header_count = header_lines
        else:
            header_count = 0
            with open(chr_file, "r", encoding="utf-8") as fh:
                while fh.readline().startswith("#"):
                    header_count += 1

        with open(chr_file, "rt", encoding="utf-8") as fh:
            if not (copy_header and i == 0):
                for _ in range(header_count):
                    fh.readline()
            shutil.copyfileobj(fh, merged_file_handler, length=1024 * 1024)

    for chr_file in file_names:
        if os.path.exists(chr_file):
            os.remove(chr_file)


def merge_counts(counter, label: str, chr_ids: Iterable[str], unaligned_reads: int = 0) -> None:
    file_name = counter.output_counts_file_name
    with counter.get_output_file_handler() as merged_file_handler:
        merge_files(file_name, label, chr_ids, merged_file_handler, header_lines=1)

        if getattr(counter, "usable_file_name", None):
            for f in merge_file_list(counter.usable_file_name, label, chr_ids):
                if os.path.exists(f):
                    counter.load_usable(f)
                    os.remove(f)

        stat_dict: dict[str, int] = defaultdict(int)
        if getattr(counter, "output_stats_file_name", None) and getattr(counter, "ignore_read_groups", True):
            stats_file_names = merge_file_list(counter.output_stats_file_name, label, chr_ids)
            for sf in stats_file_names:
                if not os.path.exists(sf):
                    continue
                with open(sf, "r", encoding="utf-8") as fh:
                    for line in fh:
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            stat_dict[parts[0]] += int(parts[1])
                os.remove(sf)

            if unaligned_reads > 0:
                stat_dict["__not_aligned"] = unaligned_reads
            for key, val in stat_dict.items():
                merged_file_handler.write(f"{key}\t{val}\n")


def normalize_path(config_path: str, file_path: str) -> str:
    if os.path.isabs(file_path):
        return os.path.normpath(file_path)
    return os.path.normpath(os.path.join(os.path.dirname(config_path), file_path))
