"""Streaming chromosome-partitioned read table splitting with strict file handle closure."""

from __future__ import annotations

import gzip
import os
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pysam

from isoquant_next.core.concurrency import WorkerPoolExecutor
from isoquant_next.core.file_naming import convert_chr_id_to_file_name_str


def get_chromosome_read_counts(bam_file: str | Path, chr_ids: Optional[Iterable[str]] = None) -> Dict[str, int]:
    """Read BAM index statistics (no full BAM scan) to determine per-chromosome mapped reads."""
    wanted = set(chr_ids) if chr_ids is not None else None
    counts: Dict[str, int] = {}
    with pysam.AlignmentFile(str(bam_file), "rb") as af:
        for stat in af.get_index_statistics():
            if wanted is None or stat.contig in wanted:
                counts[stat.contig] = int(stat.mapped)
    return counts


def split_read_table_by_chromosomes(
    table_path: str | Path,
    bam_path: str | Path,
    output_prefix: str | Path,
    chr_ids: Optional[List[str]] = None,
    has_header: bool = True,
) -> Dict[str, str]:
    """
    Split read-keyed TSV(.gz) into per-chromosome shards using bounded batch chromosome passes.
    Avoids opening thousands of simultaneous handles across 700+ scaffold contigs.
    """
    table_p = Path(table_path)
    out_base = Path(output_prefix)
    out_base.parent.mkdir(parents=True, exist_ok=True)

    with pysam.AlignmentFile(str(bam_path), "rb") as af:
        target_chroms = list(chr_ids) if chr_ids is not None else list(af.references)

    # First index read_id -> chr_id from single pass over BAM (compact bytes->int16 map)
    chrom_to_idx = {c: i for i, c in enumerate(target_chroms)}
    read_to_chr: Dict[str, int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as af:
        for aln in af.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            cid = chrom_to_idx.get(aln.reference_name or "")
            if cid is not None and aln.query_name is not None:
                read_to_chr[aln.query_name] = cid

    opener = gzip.open if table_p.suffix == ".gz" else open
    out_paths: Dict[str, str] = {}
    # Write in manageable FD batches (max 64 open files at a time)
    buffers: Dict[int, List[str]] = defaultdict(list)
    seen_header: Optional[str] = None

    def _flush_cid(cid: int, lines: List[str]) -> None:
        if not lines:
            return
        chr_name = target_chroms[cid]
        out_f = out_paths.get(chr_name)
        first = out_f is None
        if out_f is None:
            out_f = f"{out_base}_{convert_chr_id_to_file_name_str(chr_name)}.tsv"
            out_paths[chr_name] = out_f
        with open(out_f, "w" if first else "a", encoding="utf-8") as fh:
            if first and seen_header is not None:
                fh.write(seen_header)
            fh.writelines(lines)
        lines.clear()

    with opener(table_p, "rt", encoding="utf-8") as inf:
        for lineno, line in enumerate(inf):
            if lineno == 0 and has_header and (line.startswith("#") or line.startswith("read_id\t")):
                seen_header = line
                continue
            if not line or line[0] == "#":
                continue
            tab_idx = line.find("\t")
            rid = line[:tab_idx] if tab_idx != -1 else line.strip()
            cid = read_to_chr.get(rid)
            if cid is not None:
                buf = buffers[cid]
                buf.append(line)
                if len(buf) >= 4096:
                    _flush_cid(cid, buf)

    for cid, buf in buffers.items():
        if buf:
            _flush_cid(cid, buf)
    return out_paths
