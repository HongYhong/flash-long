"""Zero-pickle CompactMultimapTuple resolution manager across worker chromosomes."""

from __future__ import annotations

from collections import defaultdict
import logging
import os
from typing import Dict, Iterable, List, Optional, Tuple

from isoquant_next.assignment.loader import BasicReadAssignmentLoader
from isoquant_next.assignment.multimap import MultimapResolver, MultimapResolvingStrategy
from isoquant_next.assignment.types import BasicReadAssignment, CompactMultimapTuple, ReadAssignmentType
from isoquant_next.core.file_naming import dynamic_pools_file_name, multimappers_file_name, saves_file_name
from isoquant_next.core.serialization import TERMINATION_INT, write_int, write_list
from isoquant_next.core.string_pools import StringPoolManager

logger = logging.getLogger("IsoQuant")


def _compact_to_basic(c: CompactMultimapTuple, string_pools: StringPoolManager) -> BasicReadAssignment:
    """Hydrate CompactMultimapTuple into BasicReadAssignment right before binary disk write."""
    b = BasicReadAssignment.__new__(BasicReadAssignment)
    b.assignment_id = c.assignment_id
    b.read_id = c.read_id
    b._string_pools = string_pools
    b.chr_id_int = c.chr_id_int
    b.barcode_id = None
    b.umi_id = None
    b.read_group_ids = []
    b.start = c.start
    b.end = c.end
    b.genomic_region = (c.genomic_start, c.genomic_end)
    b.multimapper = c.multimapper
    b.polyA_found = c.polyA_found
    b.assignment_type = c.assignment_type
    b.gene_assignment_type = c.gene_assignment_type
    b.penalty_score = c.penalty_score
    b.gene_ids = list(c.gene_ids)
    b.isoform_ids = list(c.isoform_ids)
    return b


def resolve_compact_multimapper_records(
    chr_ids: Iterable[str],
    sample,
    multimapped_reads: Dict[str, List[CompactMultimapTuple]],
    strategy: MultimapResolvingStrategy,
    string_pools: StringPoolManager,
) -> Tuple[int, int]:
    """Resolve CompactMultimapTuple clusters and emit per-chr `.multimappers_<chr>` binary shards."""
    resolver = MultimapResolver(strategy)
    chr_id_list = list(chr_ids)
    dumpers = {cid: open(multimappers_file_name(sample.out_raw_file, cid), "wb") for cid in chr_id_list}
    total_assignments = 0
    polya_assignments = 0
    try:
        for comp_list in multimapped_reads.values():
            if len(comp_list) > 1:
                resolved = resolver.resolve(comp_list) or comp_list
                resolved_by_chr: Dict[str, List[BasicReadAssignment]] = defaultdict(list)
                for c in resolved:
                    chr_name = (
                        string_pools.chromosome_pool.get_str(c.chr_id_int)
                        if c.chr_id_int >= 0
                        else chr_id_list[0]
                    )
                    resolved_by_chr[chr_name].append(_compact_to_basic(c, string_pools))
                for chr_name, blist in resolved_by_chr.items():
                    if chr_name in dumpers:
                        write_list(blist, dumpers[chr_name], BasicReadAssignment.serialize)
                out_list = resolved
            else:
                out_list = comp_list

            for a in out_list:
                if a.assignment_type != ReadAssignmentType.suspended:
                    total_assignments += 1
                    if a.polyA_found:
                        polya_assignments += 1

        for cid in chr_id_list:
            write_int(TERMINATION_INT, dumpers[cid])
    finally:
        for fh in dumpers.values():
            fh.close()
    return total_assignments, polya_assignments


def prepare_compact_multimappers_from_saves(
    chr_ids: Iterable[str],
    sample,
    multimappers_counts: Dict[str, int],
    string_pools: StringPoolManager,
) -> Tuple[Dict[str, List[CompactMultimapTuple]], int, int]:
    """Scan `.save_<chr>` once for multi-hit read IDs and accumulate compact 56-byte tuples."""
    multimapped_reads: Dict[str, List[CompactMultimapTuple]] = defaultdict(list)
    unique_assignments = 0
    polya_unique_assignments = 0

    for chr_id in chr_ids:
        dyn_file = dynamic_pools_file_name(sample.out_raw_file, chr_id)
        if os.path.exists(dyn_file):
            with open(dyn_file, "rb") as fh:
                string_pools.deserialize_dynamic_pools(fh)
        chr_dump_file = saves_file_name(sample.out_raw_file, chr_id)
        with BasicReadAssignmentLoader(chr_dump_file, string_pools) as loader:
            while loader.has_next():
                for ra in loader.get_next():
                    if ra is None:
                        continue
                    if multimappers_counts.get(ra.read_id, 1) == 1:
                        unique_assignments += 1
                        if ra.polyA_found:
                            polya_unique_assignments += 1
                        continue
                    multimapped_reads[ra.read_id].append(ra.to_compact_tuple())
    return multimapped_reads, unique_assignments, polya_unique_assignments


class ProcessedReadsCollector:
    """Per-worker summary returned instead of heavyweight Python assignment object graphs."""

    __slots__ = ("total_assignments", "polya_assignments", "compact_payloads", "read_id_counts")

    def __init__(self) -> None:
        self.total_assignments = 0
        self.polya_assignments = 0
        self.compact_payloads: List[bytes] = []
        self.read_id_counts: List[Tuple[str, int]] = []
