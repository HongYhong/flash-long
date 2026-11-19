"""K-mer indexers including NEP 50 uint64-clean SharedMemoryArray2BitKmerIndexer."""

from __future__ import annotations

import logging
import math
import sys
from collections import defaultdict
from multiprocessing import shared_memory
from typing import DefaultDict, Dict, Iterable, List, Optional, Tuple

import numpy as np

from isoquant_next.barcodes.common import bit_to_str, str_to_2bit

logger = logging.getLogger("IsoQuant")


class KmerIndexer:
    """Dictionary-backed k-mer substring indexer for small reference sets."""

    def __init__(self, known_strings: Iterable[str], kmer_size: int = 6) -> None:
        self.seq_list: List[str] = list(known_strings)
        self.k = kmer_size
        self.index: DefaultDict[str, List[int]] = defaultdict(list)
        for idx, s in enumerate(self.seq_list):
            for kmer in self._get_kmers(s):
                self.index[kmer].append(idx)

    def _get_kmers(self, seq: str) -> Iterable[str]:
        if len(seq) < self.k:
            return
        for i in range(len(seq) - self.k + 1):
            yield seq[i : i + self.k]

    def _get_kmers_substr(self, seq: str, start: int, end: int) -> Iterable[str]:
        if end - start < self.k:
            return
        for i in range(start, end - self.k + 1):
            yield seq[i : i + self.k]

    def append(self, barcode: str) -> None:
        idx = len(self.seq_list)
        self.seq_list.append(barcode)
        for kmer in self._get_kmers(barcode):
            self.index[kmer].append(idx)

    def empty(self) -> bool:
        return len(self.seq_list) == 0

    def get_occurrences(
        self,
        sequence: str,
        max_hits: int = 0,
        min_kmers: int = 1,
        hits_delta: int = 1,
        ignore_equal: bool = False,
    ) -> List[Tuple[str, int, List[int]]]:
        counts: DefaultDict[int, int] = defaultdict(int)
        positions: DefaultDict[int, List[int]] = defaultdict(list)
        for pos, kmer in enumerate(self._get_kmers(sequence)):
            for i in self.index.get(kmer, ()):
                counts[i] += 1
                positions[i].append(pos)

        result = [
            (self.seq_list[i], cnt, positions[i])
            for i, cnt in counts.items()
            if cnt >= min_kmers and not (ignore_equal and self.seq_list[i] == sequence)
        ]
        if not result:
            return []
        top_hits = max(result, key=lambda x: x[1])[1]
        filtered = sorted((r for r in result if r[1] >= top_hits - hits_delta), key=lambda x: x[1], reverse=True)
        return filtered if max_hits == 0 else filtered[:max_hits]

    def get_occurrences_substr(
        self,
        sequence: str,
        start: int,
        end: int,
        max_hits: int = 0,
        min_kmers: int = 1,
        hits_delta: int = 1,
    ) -> List[Tuple[str, int, List[int]]]:
        counts: DefaultDict[int, int] = defaultdict(int)
        positions: DefaultDict[int, List[int]] = defaultdict(list)
        for pos, kmer in enumerate(self._get_kmers_substr(sequence, start, end)):
            for i in self.index.get(kmer, ()):
                counts[i] += 1
                positions[i].append(pos)
        result = [(self.seq_list[i], cnt, positions[i]) for i, cnt in counts.items() if cnt >= min_kmers]
        if not result:
            return []
        top_hits = max(result, key=lambda x: x[1])[1]
        filtered = sorted((r for r in result if r[1] >= top_hits - hits_delta), key=lambda x: x[1], reverse=True)
        return filtered if max_hits == 0 else filtered[:max_hits]


class ArrayKmerIndexer:
    """2-bit dense 4^k table indexer for k <= 8."""

    NUCL2BIN: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}

    def __init__(self, known_strings: Iterable[str], kmer_size: int = 6) -> None:
        self.seq_list: List[str] = list(known_strings)
        self.k = kmer_size
        total_kmers = 1 << (2 * kmer_size)
        self.index: List[List[int]] = [[] for _ in range(total_kmers)]
        self.mask = (1 << (2 * self.k)) - 1
        for i, bc in enumerate(self.seq_list):
            for kmer_idx in self._get_kmer_indexes(bc):
                self.index[kmer_idx].append(i)

    def _get_kmer_indexes(self, seq: str) -> Iterable[int]:
        if len(seq) < self.k:
            return
        kmer_idx = 0
        for i in range(self.k):
            kmer_idx |= ArrayKmerIndexer.NUCL2BIN[seq[i]] << ((self.k - i - 1) * 2)
        yield kmer_idx
        for i in range(self.k, len(seq)):
            kmer_idx = ((kmer_idx << 2) & self.mask) | ArrayKmerIndexer.NUCL2BIN[seq[i]]
            yield kmer_idx

    def append(self, barcode: str) -> None:
        idx = len(self.seq_list)
        self.seq_list.append(barcode)
        for kmer_idx in self._get_kmer_indexes(barcode):
            self.index[kmer_idx].append(idx)

    def empty(self) -> bool:
        return len(self.seq_list) == 0

    def get_occurrences(
        self,
        sequence: str,
        max_hits: int = 0,
        min_kmers: int = 1,
        hits_delta: int = 1,
        ignore_equal: bool = False,
    ) -> List[Tuple[str, int, List[int]]]:
        counts: DefaultDict[int, int] = defaultdict(int)
        positions: DefaultDict[int, List[int]] = defaultdict(list)
        for pos, kmer_idx in enumerate(self._get_kmer_indexes(sequence)):
            for i in self.index[kmer_idx]:
                counts[i] += 1
                positions[i].append(pos)
        result = [
            (self.seq_list[i], cnt, positions[i])
            for i, cnt in counts.items()
            if cnt >= min_kmers and not (ignore_equal and self.seq_list[i] == sequence)
        ]
        if not result:
            return []
        top_hits = max(result, key=lambda x: x[1])[1]
        filtered = sorted((r for r in result if r[1] >= top_hits - hits_delta), key=lambda x: x[1], reverse=True)
        return filtered if max_hits == 0 else filtered[:max_hits]


class Dict2BitKmerIndexer:
    """Memory-efficient dictionary k-mer index over 2-bit integer barcodes (A=0,C=1,T=2,G=3)."""

    def __init__(self, known_bit_seqs: Iterable[int], kmer_size: int, seq_len: int) -> None:
        self.seq_list: List[int] = list(known_bit_seqs)
        self.k = kmer_size
        self.seq_len = seq_len
        self.mask = (1 << (2 * self.k)) - 1
        self.index: DefaultDict[int, List[int]] = defaultdict(list)
        for i, bin_seq in enumerate(self.seq_list):
            for kmer_idx in self._get_kmer_indexes_from_bits(bin_seq):
                self.index[kmer_idx].append(i)

    def _get_kmer_indexes_from_bits(self, bin_seq: int) -> Iterable[int]:
        v = int(bin_seq)
        for i in range(self.seq_len - self.k + 1):
            yield (v >> ((self.seq_len - self.k - i) * 2)) & self.mask

    def _get_kmer_indexes_from_str(self, seq: str) -> Iterable[int]:
        if len(seq) < self.k:
            return
        kmer_idx = 0
        for i in range(self.k):
            kmer_idx |= ((ord(seq[i]) & 6) >> 1) << ((self.k - i - 1) * 2)
        yield kmer_idx
        for i in range(self.k, len(seq)):
            kmer_idx = ((kmer_idx << 2) & self.mask) | ((ord(seq[i]) & 6) >> 1)
            yield kmer_idx

    def empty(self) -> bool:
        return len(self.seq_list) == 0

    def get_occurrences(
        self,
        sequence: str,
        max_hits: int = 0,
        min_kmers: int = 1,
        hits_delta: int = 1,
        ignore_equal: bool = False,
    ) -> List[Tuple[str, int, List[int]]]:
        counts: DefaultDict[int, int] = defaultdict(int)
        positions: DefaultDict[int, List[int]] = defaultdict(list)
        query_bin = str_to_2bit(sequence)
        for pos, kmer_idx in enumerate(self._get_kmer_indexes_from_str(sequence)):
            for idx in self.index.get(kmer_idx, ()):
                counts[idx] += 1
                positions[idx].append(pos)

        result: List[Tuple[str, int, List[int]]] = []
        for idx, cnt in counts.items():
            if cnt < min_kmers:
                continue
            bc_int = int(self.seq_list[idx])
            if ignore_equal and bc_int == query_bin and len(sequence) == self.seq_len:
                continue
            result.append((bit_to_str(bc_int, self.seq_len), cnt, positions[idx]))
        if not result:
            return []
        top_hits = max(result, key=lambda x: x[1])[1]
        filtered = sorted((r for r in result if r[1] >= top_hits - hits_delta), key=lambda x: x[1], reverse=True)
        return filtered if max_hits == 0 else filtered[:max_hits]


class Array2BitKmerIndexer:
    """Flat CSR NumPy k-mer index over 2-bit integer barcodes using explicit uint64 math (NEP 50 safe)."""

    SEQ_DTYPE = np.uint64
    KMER_DTYPE = np.uint32
    INDEX_DTYPE = np.int64

    def __init__(self, known_bit_seqs, kmer_size: int, seq_len: int) -> None:
        self.k = int(kmer_size)
        self.seq_len = int(seq_len)
        self.seq_list = np.ascontiguousarray(known_bit_seqs, dtype=self.SEQ_DTYPE)
        self.mask = np.uint64((1 << (2 * self.k)) - 1)
        self.total_kmers = 1 << (2 * self.k)
        self.num_seqs = int(len(self.seq_list))

        if self.num_seqs == 0 or self.seq_len < self.k:
            self.starts = np.zeros(self.total_kmers + 1, dtype=self.INDEX_DTYPE)
            self.seq_indices = np.zeros(0, dtype=self.INDEX_DTYPE)
            return

        kmers_per_seq = self.seq_len - self.k + 1
        shifts = np.arange(kmers_per_seq - 1, -1, -1, dtype=np.uint64) * np.uint64(2)
        all_kmers = ((self.seq_list[:, None] >> shifts[None, :]) & self.mask).astype(self.KMER_DTYPE).ravel()
        seq_ids = np.repeat(np.arange(self.num_seqs, dtype=self.INDEX_DTYPE), kmers_per_seq)

        order = np.argsort(all_kmers, kind="mergesort")
        sorted_kmers = all_kmers[order]
        self.seq_indices = np.ascontiguousarray(seq_ids[order], dtype=self.INDEX_DTYPE)
        counts = np.bincount(sorted_kmers, minlength=self.total_kmers).astype(self.INDEX_DTYPE)
        self.starts = np.empty(self.total_kmers + 1, dtype=self.INDEX_DTYPE)
        self.starts[0] = 0
        np.cumsum(counts, out=self.starts[1:])

    @property
    def total_sequences(self) -> int:
        return self.num_seqs

    def _get_kmer_bin_indexes(self, seq) -> Iterable[int]:
        if isinstance(seq, (int, np.integer)):
            v = int(seq)
            if self.seq_len < self.k:
                return
            mask_py = (1 << (2 * self.k)) - 1
            for i in range(self.seq_len - self.k + 1):
                yield (v >> ((self.seq_len - self.k - i) * 2)) & mask_py
            return
        if len(seq) < self.k:
            return
        mask_py = (1 << (2 * self.k)) - 1
        kmer_idx = 0
        for i in range(self.k):
            kmer_idx |= ((ord(seq[i]) & 6) >> 1) << ((self.k - i - 1) * 2)
        yield kmer_idx
        for i in range(self.k, len(seq)):
            kmer_idx = ((kmer_idx << 2) & mask_py) | ((ord(seq[i]) & 6) >> 1)
            yield kmer_idx

    def empty(self) -> bool:
        return self.num_seqs == 0

    def get_occurrences(
        self,
        sequence: str,
        max_hits: int = 0,
        min_kmers: int = 1,
        hits_delta: int = 1,
        ignore_equal: bool = False,
    ) -> List[Tuple[str, int, List[int]]]:
        if self.num_seqs == 0:
            return []
        counts: DefaultDict[int, int] = defaultdict(int)
        positions: DefaultDict[int, List[int]] = defaultdict(list)
        query_bin = str_to_2bit(sequence) if ignore_equal else 0

        for pos, kmer_idx in enumerate(self._get_kmer_bin_indexes(sequence)):
            s0 = int(self.starts[kmer_idx])
            s1 = int(self.starts[kmer_idx + 1])
            for raw_id in self.seq_indices[s0:s1]:
                idx = int(raw_id)
                counts[idx] += 1
                positions[idx].append(pos)

        result: List[Tuple[str, int, List[int]]] = []
        for idx, cnt in counts.items():
            if cnt < min_kmers:
                continue
            bc_int = int(self.seq_list[idx])
            if ignore_equal and bc_int == query_bin and len(sequence) == self.seq_len:
                continue
            result.append((bit_to_str(bc_int, self.seq_len), cnt, positions[idx]))
        if not result:
            return []
        top_hits = max(result, key=lambda x: x[1])[1]
        filtered = sorted((r for r in result if r[1] >= top_hits - hits_delta), key=lambda x: x[1], reverse=True)
        return filtered if max_hits == 0 else filtered[:max_hits]


def _open_shm_readonly(name: str) -> shared_memory.SharedMemory:
    """Attach existing shared memory safely without resource_tracker double-unlink on Python 3.13+."""
    if sys.version_info >= (3, 13):
        return shared_memory.SharedMemory(name=name, create=False, track=False)
    return shared_memory.SharedMemory(name=name, create=False)


class SharedMemoryIndexInfo:
    """Picklable descriptor for attaching to SharedMemoryArray2BitKmerIndexer across processes."""

    def __init__(
        self,
        barcode_count: int,
        kmer_size: int,
        seq_len: int,
        index_size: int,
        barcodes_sm_name: str,
        index_sm_name: str,
        index_range_sm_name: str,
    ) -> None:
        self.barcode_count = int(barcode_count)
        self.kmer_size = int(kmer_size)
        self.seq_len = int(seq_len)
        self.index_size = int(index_size)
        self.barcodes_sm_name = barcodes_sm_name
        self.index_sm_name = index_sm_name
        self.index_range_sm_name = index_range_sm_name

    def __getstate__(self) -> tuple:
        return (
            self.barcode_count,
            self.kmer_size,
            self.seq_len,
            self.index_size,
            self.barcodes_sm_name,
            self.index_sm_name,
            self.index_range_sm_name,
        )

    def __setstate__(self, state: tuple) -> None:
        (
            self.barcode_count,
            self.kmer_size,
            self.seq_len,
            self.index_size,
            self.barcodes_sm_name,
            self.index_sm_name,
            self.index_range_sm_name,
        ) = state


class SharedMemoryArray2BitKmerIndexer:
    """Zero-copy POSIX shared-memory Array2Bit indexer with deterministic context-managed cleanup."""

    SEQ_DTYPE = np.uint64
    KMER_DTYPE = np.uint32
    INDEX_DTYPE = np.int64

    def __init__(
        self,
        known_bit_seqs=None,
        kmer_size: Optional[int] = None,
        seq_len: Optional[int] = None,
        _from_shared=None,
    ) -> None:
        self._shm_handles: List[shared_memory.SharedMemory] = []
        self._owns_shm: bool = False

        if _from_shared is not None:
            if isinstance(_from_shared, SharedMemoryIndexInfo):
                info = _from_shared
                self.k = int(info.kmer_size)
                self.seq_len = int(info.seq_len)
                self.num_seqs = int(info.barcode_count)
                self.total_kmers = 1 << (2 * self.k)
                num_entries = int(info.index_size)
                seqs_shm_name = info.barcodes_sm_name
                starts_shm_name = info.index_range_sm_name
                seq_indices_shm_name = info.index_sm_name
            else:
                self.k = int(_from_shared["k"])
                self.seq_len = int(_from_shared["seq_len"])
                self.num_seqs = int(_from_shared["num_seqs"])
                self.total_kmers = int(_from_shared["total_kmers"])
                num_entries = int(_from_shared["num_entries"])
                seqs_shm_name = _from_shared["seqs_shm"]
                starts_shm_name = _from_shared["starts_shm"]
                seq_indices_shm_name = _from_shared["seq_indices_shm"]

            self.mask = np.uint64((1 << (2 * self.k)) - 1)
            self._shm_info = SharedMemoryIndexInfo(
                barcode_count=self.num_seqs,
                kmer_size=self.k,
                seq_len=self.seq_len,
                index_size=num_entries,
                barcodes_sm_name=seqs_shm_name,
                index_sm_name=seq_indices_shm_name,
                index_range_sm_name=starts_shm_name,
            )

            shm_s = _open_shm_readonly(seqs_shm_name)
            shm_st = _open_shm_readonly(starts_shm_name)
            shm_si = _open_shm_readonly(seq_indices_shm_name)
            self._shm_handles.extend([shm_s, shm_st, shm_si])
            self.seq_list = np.ndarray((self.num_seqs,), dtype=self.SEQ_DTYPE, buffer=shm_s.buf)
            self.starts = np.ndarray((self.total_kmers + 1,), dtype=self.INDEX_DTYPE, buffer=shm_st.buf)
            self.seq_indices = np.ndarray((num_entries,), dtype=self.INDEX_DTYPE, buffer=shm_si.buf)
            return

        assert kmer_size is not None and seq_len is not None and known_bit_seqs is not None
        self._owns_shm = True
        local = Array2BitKmerIndexer(known_bit_seqs, kmer_size, seq_len)
        self.k = local.k
        self.seq_len = local.seq_len
        self.num_seqs = local.num_seqs
        self.total_kmers = local.total_kmers
        self.mask = local.mask

        def _alloc(arr: np.ndarray) -> Tuple[shared_memory.SharedMemory, np.ndarray]:
            nbytes = max(1, int(arr.nbytes))
            shm = shared_memory.SharedMemory(create=True, size=nbytes)
            self._shm_handles.append(shm)
            view = np.ndarray(arr.shape, dtype=arr.dtype, buffer=shm.buf)
            if arr.size > 0:
                view[:] = arr
            return shm, view

        shm_s, self.seq_list = _alloc(local.seq_list)
        shm_st, self.starts = _alloc(local.starts)
        shm_si, self.seq_indices = _alloc(local.seq_indices)
        self._shm_info = SharedMemoryIndexInfo(
            barcode_count=self.num_seqs,
            kmer_size=self.k,
            seq_len=self.seq_len,
            index_size=int(local.seq_indices.size),
            barcodes_sm_name=shm_s.name,
            index_sm_name=shm_si.name,
            index_range_sm_name=shm_st.name,
        )

    @property
    def total_sequences(self) -> int:
        return self.num_seqs

    @property
    def index_size(self) -> int:
        return int(self.seq_indices.size)

    def __enter__(self) -> "SharedMemoryArray2BitKmerIndexer":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.cleanup()

    @classmethod
    def from_sharable_info(cls, info) -> "SharedMemoryArray2BitKmerIndexer":
        return cls(_from_shared=info)

    def get_sharable_info(self) -> SharedMemoryIndexInfo:
        state = self._shm_info.__getstate__()
        cloned = SharedMemoryIndexInfo.__new__(SharedMemoryIndexInfo)
        cloned.__setstate__(state)
        return cloned

    def empty(self) -> bool:
        return self.num_seqs == 0

    def get_occurrences(
        self,
        sequence: str,
        max_hits: int = 0,
        min_kmers: int = 1,
        hits_delta: int = 1,
        ignore_equal: bool = False,
    ) -> List[Tuple[str, int, List[int]]]:
        if self.num_seqs == 0:
            return []
        counts: DefaultDict[int, int] = defaultdict(int)
        positions: DefaultDict[int, List[int]] = defaultdict(list)
        query_bin = str_to_2bit(sequence) if ignore_equal else 0

        mask_py = (1 << (2 * self.k)) - 1
        if len(sequence) >= self.k:
            kmer_idx = 0
            for i in range(self.k):
                kmer_idx |= ((ord(sequence[i]) & 6) >> 1) << ((self.k - i - 1) * 2)
            for pos in range(len(sequence) - self.k + 1):
                if pos > 0:
                    kmer_idx = ((kmer_idx << 2) & mask_py) | ((ord(sequence[pos + self.k - 1]) & 6) >> 1)
                s0 = int(self.starts[kmer_idx])
                s1 = int(self.starts[kmer_idx + 1])
                for raw_id in self.seq_indices[s0:s1]:
                    idx = int(raw_id)
                    counts[idx] += 1
                    positions[idx].append(pos)

        result: List[Tuple[str, int, List[int]]] = []
        for idx, cnt in counts.items():
            if cnt < min_kmers:
                continue
            bc_int = int(self.seq_list[idx])
            if ignore_equal and bc_int == query_bin and len(sequence) == self.seq_len:
                continue
            result.append((bit_to_str(bc_int, self.seq_len), cnt, positions[idx]))
        if not result:
            return []
        top_hits = max(result, key=lambda x: x[1])[1]
        filtered = sorted((r for r in result if r[1] >= top_hits - hits_delta), key=lambda x: x[1], reverse=True)
        return filtered if max_hits == 0 else filtered[:max_hits]

    def close(self) -> None:
        self.seq_list = np.empty(0, dtype=self.SEQ_DTYPE)
        self.starts = np.empty(0, dtype=self.INDEX_DTYPE)
        self.seq_indices = np.empty(0, dtype=self.INDEX_DTYPE)
        for shm in self._shm_handles:
            try:
                shm.close()
            except Exception:
                pass
        if self._owns_shm:
            for shm in self._shm_handles:
                try:
                    shm.unlink()
                except FileNotFoundError:
                    pass
                except Exception:
                    pass
        self._shm_handles.clear()

    def cleanup(self) -> None:
        self.close()
