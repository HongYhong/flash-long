"""Worker-safe lazy reference genome wrapper powered by pysam.FastaFile."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Optional

import pysam


class ReferenceGenome:
    """Lazy context-managed FASTA accessor supporting plain FASTA and bgzipped .fa.gz."""

    def __init__(self, fasta_path: str | Path, index_dir: Optional[str | Path] = None) -> None:
        self.fasta_path = Path(fasta_path)
        self.index_dir = Path(index_dir) if index_dir else None
        self._handle: Optional[pysam.FastaFile] = None
        self._lengths: Optional[Dict[str, int]] = None

    def _ensure_open(self) -> pysam.FastaFile:
        if self._handle is not None:
            return self._handle
        fai_path = Path(str(self.fasta_path) + ".fai")
        if fai_path.exists():
            self._handle = pysam.FastaFile(str(self.fasta_path), filepath_index=str(fai_path))
        else:
            # If parent directory is read-only, copy/symlink FASTA to index_dir before building .fai
            try:
                self._handle = pysam.FastaFile(str(self.fasta_path))
            except OSError:
                if self.index_dir is None:
                    raise
                self.index_dir.mkdir(parents=True, exist_ok=True)
                linked = self.index_dir / self.fasta_path.name
                if not linked.exists():
                    os.symlink(self.fasta_path.resolve(), linked)
                self._handle = pysam.FastaFile(str(linked))
        return self._handle

    def __enter__(self) -> "ReferenceGenome":
        self._ensure_open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def get_contig_lengths(self) -> Dict[str, int]:
        if self._lengths is None:
            fa = self._ensure_open()
            self._lengths = dict(zip(fa.references, fa.lengths))
        return self._lengths

    def fetch_1based(self, chrom: str, start_1based: int, end_1based: int) -> str:
        """Fetch inclusive 1-based [start, end] reference slice via 0-based half-open pysam query."""
        fa = self._ensure_open()
        s0 = max(0, start_1based - 1)
        e0 = max(s0, end_1based)
        return fa.fetch(chrom, s0, e0)

    def fetch_0based(self, chrom: str, start_0based: int, end_0based: int) -> str:
        fa = self._ensure_open()
        return fa.fetch(chrom, max(0, start_0based), max(start_0based, end_0based))

    def contig_view(self, chrom: str) -> "LazyContigView":
        """Return a zero-memory sliceable proxy for ``chrom`` backed by pysam.FastaFile."""
        return LazyContigView(self, chrom, self.get_contig_lengths().get(chrom, 0))


class LazyContigView:
    """Sliceable O(1)-memory chromosome proxy avoiding 250 MiB full-chromosome string loads."""

    __slots__ = ("_ref", "_chrom", "_length")

    def __init__(self, ref: ReferenceGenome, chrom: str, length: int) -> None:
        self._ref = ref
        self._chrom = chrom
        self._length = int(length)

    def __len__(self) -> int:
        return self._length

    def __getitem__(self, key):
        if isinstance(key, slice):
            start, stop, step = key.indices(self._length)
            if step != 1:
                return self._ref.fetch_0based(self._chrom, 0, self._length)[slice(start, stop, step)]
            if stop <= start:
                return ""
            return self._ref.fetch_0based(self._chrom, start, stop)
        if key < 0:
            key += self._length
        if key < 0 or key >= self._length:
            raise IndexError("LazyContigView index out of range")
        return self._ref.fetch_0based(self._chrom, key, key + 1)
