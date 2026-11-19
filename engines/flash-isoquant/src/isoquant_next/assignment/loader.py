"""Resource-safe binary assignment loaders backed by `.iqidx` mmap and LazyContigView."""

from __future__ import annotations

from collections import defaultdict
import logging
import os
from typing import Dict, Optional, Set

from isoquant_next.annotation.gene_info import GeneInfo
from isoquant_next.annotation.iqidx import MmapAnnotationIndex
from isoquant_next.annotation.reference import ReferenceGenome
from isoquant_next.assignment.io import (
    GeneListTmpFileAssignmentLoader,
    NormalTmpFileAssignmentLoader,
    QuickTmpFileAssignmentLoader,
)
from isoquant_next.assignment.types import BasicReadAssignment, ReadAssignmentType
from isoquant_next.core.file_naming import (
    filtered_reads_file_name,
    multimappers_file_name,
    saves_file_name,
)
from isoquant_next.core.serialization import TERMINATION_INT, read_int

logger = logging.getLogger("IsoQuant")


def open_genedb_index(genedb):
    """Open `.iqidx` via zero-copy mmap or pass through an already-opened index/DB."""
    if genedb is None:
        return None
    if isinstance(genedb, str):
        if genedb.endswith(".iqidx"):
            return MmapAnnotationIndex(genedb)
        raise ValueError(f"Expected .iqidx annotation index path, got: {genedb}")
    return genedb


class BasicAssignmentLoader:
    """Base binary assignment stream loader."""

    def __init__(self, save_file_name: str) -> None:
        logger.info("Loading read assignments from %s", save_file_name)
        assert os.path.exists(save_file_name), f"Missing save file: {save_file_name}"
        self.save_file_name = save_file_name

    def close(self) -> None:
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def has_next(self) -> bool:
        raise NotImplementedError()

    def get_next(self):
        raise NotImplementedError()


class FullAssignmentLoader(BasicAssignmentLoader):
    """Applies multimap resolutions and optional UMI read filter sets on stream load."""

    def __init__(self, save_file_name: str, multimapped_chr_dict, filtered_read_set=None) -> None:
        super().__init__(save_file_name)
        self.multimapped_chr_dict = multimapped_chr_dict
        self.filtered_read_set = filtered_read_set

    def load_read_assignment(self, read_assignment):
        if self.filtered_read_set is not None and read_assignment.read_id not in self.filtered_read_set:
            return None
        if self.multimapped_chr_dict is not None and read_assignment.read_id in self.multimapped_chr_dict:
            resolved_assignment = None
            for cand in self.multimapped_chr_dict[read_assignment.read_id]:
                if cand.assignment_id == read_assignment.assignment_id and (
                    getattr(cand, "chr_id", read_assignment.chr_id) == read_assignment.chr_id
                ):
                    resolved_assignment = cand
            if resolved_assignment is None:
                return None
            if resolved_assignment.assignment_type == ReadAssignmentType.suspended:
                return None
            read_assignment.assignment_type = resolved_assignment.assignment_type
            read_assignment.gene_assignment_type = resolved_assignment.gene_assignment_type
            read_assignment.multimapper = resolved_assignment.multimapper
        return read_assignment


class ReadAssignmentLoader(FullAssignmentLoader):
    """Streams GeneInfo blocks + ReadAssignment lists using LazyContigView and MmapAnnotationIndex."""

    def __init__(
        self,
        save_file_name: str,
        genedb,
        chr_record,
        multimapped_chr_dict,
        string_pools,
        filtered_read_set=None,
        *,
        _owned_db=None,
        _owned_ref: Optional[ReferenceGenome] = None,
    ) -> None:
        super().__init__(save_file_name, multimapped_chr_dict, filtered_read_set)
        self.genedb = genedb
        self.chr_record = chr_record
        self._owned_db = _owned_db
        self._owned_ref = _owned_ref
        self.unpickler = NormalTmpFileAssignmentLoader(save_file_name, genedb, chr_record, string_pools)

    def close(self) -> None:
        if getattr(self, "unpickler", None) is not None:
            self.unpickler.close()
        if getattr(self, "_owned_db", None) is not None and hasattr(self._owned_db, "close"):
            self._owned_db.close()
            self._owned_db = None
        if getattr(self, "_owned_ref", None) is not None:
            self._owned_ref.close()
            self._owned_ref = None

    def has_next(self) -> bool:
        return self.unpickler.has_next()

    def get_next(self):
        if not self.unpickler.has_next():
            return None, None
        assert self.unpickler.is_gene_info()
        gene_info = self.unpickler.get_object()
        assignment_storage = []
        while self.unpickler.is_read_assignment():
            read_assignment = self.load_read_assignment(self.unpickler.get_object())
            if read_assignment is not None:
                assignment_storage.append(read_assignment)
        return gene_info, assignment_storage


class MergingSimpleReadAssignmentLoader(FullAssignmentLoader):
    """Merges adjacent overlapping gene clusters without materializing full GeneInfo."""

    def __init__(self, save_file_name: str, multimapped_chr_dict, string_pools, filtered_read_set=None) -> None:
        super().__init__(save_file_name, multimapped_chr_dict, filtered_read_set)
        self.unpickler = GeneListTmpFileAssignmentLoader(save_file_name, string_pools)
        self.current_gene_list = None

    def close(self) -> None:
        if getattr(self, "unpickler", None) is not None:
            self.unpickler.close()

    def has_next(self) -> bool:
        return self.unpickler.has_next()

    def get_next(self):
        if not self.unpickler.has_next():
            return None, None
        assignment_storage = []
        while self.unpickler.has_next():
            if self.current_gene_list is None:
                assert self.unpickler.is_gene_info()
                self.current_gene_list = self.unpickler.get_object()
            elif self.unpickler.is_gene_info():
                gene_list = self.unpickler.get_object()
                if self.current_gene_list.overlaps(gene_list):
                    self.current_gene_list.merge(gene_list)
                else:
                    self.current_gene_list = gene_list
                    return None, assignment_storage

            while self.unpickler.is_read_assignment():
                read_assignment = self.load_read_assignment(self.unpickler.get_object())
                if read_assignment is not None:
                    assignment_storage.append(read_assignment)
        return None, assignment_storage


class BasicReadAssignmentLoader(BasicAssignmentLoader):
    """Streams lightweight BasicReadAssignment records for multimap scanning."""

    def __init__(self, save_file_name: str, string_pools) -> None:
        super().__init__(save_file_name)
        self.unpickler = QuickTmpFileAssignmentLoader(save_file_name, string_pools)

    def close(self) -> None:
        if getattr(self, "unpickler", None) is not None:
            self.unpickler.close()

    def has_next(self) -> bool:
        return self.unpickler.has_next()

    def get_next(self):
        if not self.unpickler.has_next():
            return
        assert self.unpickler.is_gene_info()
        self.unpickler.get_object()
        while self.unpickler.is_read_assignment():
            yield self.unpickler.get_object()


def prepare_multimapped_reads(saves_prefix: str, chr_id: str, string_pools):
    multimapped_reads = defaultdict(list)
    mm_path = multimappers_file_name(saves_prefix, chr_id)
    if not os.path.exists(mm_path):
        return multimapped_reads
    with open(mm_path, "rb") as multimap_loader:
        list_size = read_int(multimap_loader)
        while list_size != TERMINATION_INT:
            for _ in range(list_size):
                a = BasicReadAssignment.deserialize(multimap_loader, string_pools)
                if a.chr_id == chr_id:
                    multimapped_reads[a.read_id].append(a)
            list_size = read_int(multimap_loader)
    return multimapped_reads


def prepare_read_filter(chr_id: str, saves_prefix: str, use_filtered_reads: bool) -> Optional[Set[str]]:
    if not use_filtered_reads:
        return None
    fpath = filtered_reads_file_name(saves_prefix, chr_id)
    if not os.path.exists(fpath):
        return None
    with open(fpath, "r", encoding="utf-8") as fh:
        return {line.rstrip() for line in fh}


def create_assignment_loader(
    chr_id: str,
    saves_prefix: str,
    genedb,
    reference_fasta: str,
    reference_fai: Optional[str],
    string_pools,
    use_filtered_reads: bool = False,
) -> ReadAssignmentLoader:
    ref_dir = os.path.dirname(reference_fai) if reference_fai else None
    ref = ReferenceGenome(reference_fasta, index_dir=ref_dir)
    current_chr_record = ref.contig_view(chr_id)
    multimapped_reads = prepare_multimapped_reads(saves_prefix, chr_id, string_pools)
    filtered_reads = prepare_read_filter(chr_id, saves_prefix, use_filtered_reads)
    owned_db = None
    db_obj = genedb
    if isinstance(genedb, str):
        db_obj = open_genedb_index(genedb)
        owned_db = db_obj
    chr_dump_file = saves_file_name(saves_prefix, chr_id)
    return ReadAssignmentLoader(
        chr_dump_file,
        db_obj,
        current_chr_record,
        multimapped_reads,
        string_pools,
        filtered_reads,
        _owned_db=owned_db,
        _owned_ref=ref,
    )


def create_merging_assignment_loader(
    chr_id: str, saves_prefix: str, string_pools, use_filtered_reads: bool = False
) -> MergingSimpleReadAssignmentLoader:
    multimapped_reads = prepare_multimapped_reads(saves_prefix, chr_id, string_pools)
    filtered_reads = prepare_read_filter(chr_id, saves_prefix, use_filtered_reads)
    chr_dump_file = saves_file_name(saves_prefix, chr_id)
    return MergingSimpleReadAssignmentLoader(chr_dump_file, multimapped_reads, string_pools, filtered_reads)
