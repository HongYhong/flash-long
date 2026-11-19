"""Deterministic string-to-integer pools replacing repeated high-cardinality Python strings."""

from __future__ import annotations

import logging
import os
from typing import BinaryIO, Dict, Iterable, List, Optional, Union

from isoquant_next.core.serialization import read_int, read_string, write_int, write_string

logger = logging.getLogger("IsoQuant")
DEFAULT_GROUP_ID = "NA"


class StringPool:
    """Bi-directional deterministic string <-> dense 0-based integer registry."""

    def __init__(self) -> None:
        self._str_to_id: Dict[str, int] = {}
        self._id_to_str: List[str] = []

    @property
    def next_id(self) -> int:
        return len(self._id_to_str)

    def add(self, value: str) -> int:
        existing = self._str_to_id.get(value)
        if existing is not None:
            return existing
        new_id = len(self._id_to_str)
        self._str_to_id[value] = new_id
        self._id_to_str.append(value)
        return new_id

    def get_int(self, value: str) -> int:
        return self.add(value)

    def get_str(self, idx: int) -> str:
        if idx < 0 or idx >= len(self._id_to_str):
            raise IndexError(f"StringPool index {idx} out of range [0, {len(self._id_to_str)})")
        return self._id_to_str[idx]

    def __len__(self) -> int:
        return len(self._id_to_str)

    def __contains__(self, value: object) -> bool:
        return value in self._str_to_id

    @property
    def str_to_int(self) -> Dict[str, int]:
        return self._str_to_id

    def __iter__(self):
        return iter(self._id_to_str)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, StringPool):
            return False
        return self._id_to_str == other._id_to_str

    def serialize(self, outf: BinaryIO) -> None:
        write_int(len(self._id_to_str), outf)
        for s in self._id_to_str:
            write_string(s, outf)

    @classmethod
    def deserialize(cls, inf: BinaryIO) -> "StringPool":
        pool = cls()
        n = read_int(inf)
        for _ in range(n):
            pool.add(read_string(inf))
        return pool


class StringPoolManager:
    """Manages annotation, barcode, UMI, chromosome, and read-group pools deterministically."""

    _GFFUTILS_TX_TYPES = (
        "transcript", "RNA", "mRNA", "miRNA", "ncRNA", "tRNA", "rRNA",
        "snRNA", "lnc_RNA", "snoRNA", "telomerase_RNA", "Y_RNA", "scRNA",
    )

    def __init__(self) -> None:
        self.gene_pool = StringPool()
        self.transcript_pool = StringPool()
        self.chromosome_pool = StringPool()
        self.barcode_pool = StringPool()
        self.umi_pool = StringPool()
        self.read_group_tsv_pools: Dict[str, StringPool] = {}
        self.read_group_dynamic_pools: Dict[int, StringPool] = {}
        self.file_name_pool = StringPool()
        self.barcode_spot_pools: Dict[int, StringPool] = {}
        self.group_spec_pool_types: Dict[int, str] = {}

    def build_from_iqidx(self, iqidx) -> None:
        """Populate gene, transcript, and chromosome pools directly from mmap annotation index."""
        for gid in sorted(iqidx.iter_gene_ids()):
            self.gene_pool.add(gid)
        for tid in sorted(iqidx.iter_transcript_ids()):
            self.transcript_pool.add(tid)
        for cid in sorted(iqidx.chrom_names):
            self.chromosome_pool.add(cid)

    def build_from_gffutils(self, gffutils_db) -> None:
        gene_ids = {g.id for g in gffutils_db.features_of_type("gene")}
        for gid in sorted(gene_ids):
            self.gene_pool.add(gid)

        tx_ids = set()
        for ftype in self._GFFUTILS_TX_TYPES:
            try:
                for tx in gffutils_db.features_of_type(ftype):
                    tx_ids.add(tx.id)
            except Exception:
                pass
        for tid in sorted(tx_ids):
            self.transcript_pool.add(tid)

    def build_from_gene_db(self, gene_db) -> None:
        for gid in sorted(gene_db.genes.keys()):
            self.gene_pool.add(gid)
        for tid in sorted(gene_db.transcripts.keys()):
            self.transcript_pool.add(tid)

    def build_chromosome_pool(self, chr_ids: Iterable[str]) -> None:
        for chr_id in sorted(chr_ids):
            self.chromosome_pool.add(chr_id)

    def build_file_name_pool(self, sample) -> None:
        file_names = set()
        if getattr(sample, "readable_names_dict", None):
            file_names.update(sample.readable_names_dict.values())
        else:
            for lib in sample.file_list:
                file_names.add(os.path.splitext(os.path.basename(lib[0]))[0])
        for name in sorted(file_names):
            self.file_name_pool.add(name)

    def build_barcode_spot_pool(
        self, barcode2spot_spec: str, column_index: int = 0, file_column: int = 1
    ) -> None:
        bc2spot_file = barcode2spot_spec.split(":")[0]
        spots = {DEFAULT_GROUP_ID}
        with open(bc2spot_file, "r", encoding="utf-8") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split("\t")
                if len(parts) > file_column:
                    spots.add(parts[file_column])
        pool = StringPool()
        for spot in sorted(spots):
            pool.add(spot)
        self.barcode_spot_pools[column_index] = pool

    def load_barcode_pool(self, barcode_file: str) -> None:
        if not os.path.exists(barcode_file):
            return
        barcodes = {DEFAULT_GROUP_ID}
        umis = set()
        with open(barcode_file, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("read_id\t"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) >= 2 and parts[1]:
                    barcodes.add(parts[1])
                if len(parts) >= 3 and parts[2]:
                    umis.add(parts[2])
        for bc in sorted(barcodes):
            self.barcode_pool.add(bc)
        for umi in sorted(umis):
            self.umi_pool.add(umi)

    def load_read_group_tsv_pool(
        self, read_group_file: str, pool_key: str, col_index: int = 1, delimiter: str = "\t"
    ) -> None:
        if not os.path.exists(read_group_file):
            return
        values = {DEFAULT_GROUP_ID}
        with open(read_group_file, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                s = line.strip()
                if not s:
                    continue
                parts = s.split(delimiter)
                if len(parts) > col_index:
                    values.add(parts[col_index])
        pool = StringPool()
        for v in sorted(values):
            pool.add(v)
        self.read_group_tsv_pools[pool_key] = pool

    def set_group_spec_pool_type(self, spec_index: int, pool_type: str) -> None:
        self.group_spec_pool_types[spec_index] = pool_type

    def _get_or_create_dynamic_pool(self, spec_index: int) -> StringPool:
        pool = self.read_group_dynamic_pools.get(spec_index)
        if pool is None:
            pool = StringPool()
            self.read_group_dynamic_pools[spec_index] = pool
        return pool

    def get_read_group_pool(self, spec_index: int) -> StringPool:
        pool_type = self.group_spec_pool_types.get(spec_index)
        if pool_type is None:
            return self._get_or_create_dynamic_pool(spec_index)
        if pool_type == "file_name":
            return self.file_name_pool
        if pool_type.startswith("barcode_spot:"):
            col_idx = int(pool_type.split(":")[1])
            return self.barcode_spot_pools.get(col_idx) or self._get_or_create_dynamic_pool(spec_index)
        if pool_type == "barcode_spot":
            return self.barcode_spot_pools.get(0) or self._get_or_create_dynamic_pool(spec_index)
        if pool_type == "barcode":
            return self.barcode_pool
        if pool_type.startswith("tsv:"):
            parts = pool_type.split(":")
            pool_key = f"{parts[1]}:{parts[2]}" if len(parts) >= 3 else parts[1]
            return self.read_group_tsv_pools.get(pool_key) or self._get_or_create_dynamic_pool(spec_index)
        return self._get_or_create_dynamic_pool(spec_index)

    def read_group_to_ids(self, group_strings: List[Optional[str]]) -> List[int]:
        if not group_strings:
            return []
        ids: List[int] = []
        for idx, gs in enumerate(group_strings):
            if gs is not None:
                ids.append(self.get_read_group_pool(idx).get_int(gs))
            else:
                ids.append(-1)
        return ids

    def read_group_from_ids(self, group_ids: List[int]) -> List[Optional[str]]:
        if not group_ids:
            return []
        out: List[Optional[str]] = []
        for idx, gid in enumerate(group_ids):
            if gid < 0:
                out.append(None)
            else:
                pool = self.get_read_group_pool(idx)
                try:
                    out.append(pool.get_str(gid))
                except IndexError:
                    out.append(str(gid))
        return out

    def has_dynamic_pools(self) -> bool:
        return any(len(p) > 0 for p in self.read_group_dynamic_pools.values())

    def serialize_dynamic_pools(self, outfile: BinaryIO) -> None:
        non_empty = {idx: p for idx, p in self.read_group_dynamic_pools.items() if len(p) > 0}
        write_int(len(non_empty), outfile)
        for spec_idx in sorted(non_empty):
            write_int(spec_idx, outfile)
            non_empty[spec_idx].serialize(outfile)

    def deserialize_dynamic_pools(self, infile: BinaryIO) -> None:
        num_specs = read_int(infile)
        for _ in range(num_specs):
            spec_idx = read_int(infile)
            self.read_group_dynamic_pools[spec_idx] = StringPool.deserialize(infile)

    def save_barcode_umi_pools(self, filepath: str) -> None:
        with open(filepath, "wb") as fh:
            self.barcode_pool.serialize(fh)
            self.umi_pool.serialize(fh)

    def load_barcode_umi_pools(self, filepath: str) -> None:
        if not os.path.exists(filepath):
            return
        with open(filepath, "rb") as fh:
            loaded_bc = StringPool.deserialize(fh)
            loaded_umi = StringPool.deserialize(fh)
        for s in loaded_bc:
            self.barcode_pool.add(s)
        for s in loaded_umi:
            self.umi_pool.add(s)

    def get_stats(self) -> str:
        lines = [
            "String Pool Statistics:",
            f"  Genes: {len(self.gene_pool)}",
            f"  Transcripts: {len(self.transcript_pool)}",
            f"  Chromosomes: {len(self.chromosome_pool)}",
            f"  Barcodes: {len(self.barcode_pool)}",
            f"  UMIs: {len(self.umi_pool)}",
            f"  File names: {len(self.file_name_pool)}",
            f"  Barcode spots: {sum(len(p) for p in self.barcode_spot_pools.values())} "
            f"(across {len(self.barcode_spot_pools)} pools)",
        ]
        for idx, pool in self.read_group_dynamic_pools.items():
            lines.append(f"  Read group dynamic (spec {idx}): {len(pool)}")
        for key, pool in self.read_group_tsv_pools.items():
            lines.append(f"  Read group TSV ({key}): {len(pool)}")
        return "\n".join(lines)
