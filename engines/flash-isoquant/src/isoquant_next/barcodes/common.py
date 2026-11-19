"""Barcode sequence utilities, SSW-equivalent alignment, and FastxFile streaming readers."""

from __future__ import annotations

import copy
import gzip
import logging
import math
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import numpy as np
import pysam

from isoquant_next.barcodes._cy_barcode_umi import (
    encode_2bit_u8_to_u64,
    find_polyt_start_bytes,
    levenshtein_bytes,
    sw_align_bytes,
)
from isoquant_next.core.error_codes import IsoQuantExitCode
from isoquant_next.core.utils import reverse_complement

logger = logging.getLogger("IsoQuant")

NUCL2BIN: Dict[str, int] = {"A": 0, "C": 1, "G": 3, "T": 2, "a": 0, "c": 1, "g": 3, "t": 2}
BIN2NUCL: List[str] = ["A", "C", "T", "G"]
base_comp: Dict[str, str] = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", " ": " "}


def reverese_complement(my_seq: str) -> str:
    """Legacy spelling alias delegating to unified IUPAC str.translate reverse_complement."""
    if " " in my_seq:
        my_seq_trans = my_seq.translate(str.maketrans("ACGTNacgtn ", "TGCANtgcan "))
        return my_seq_trans[::-1]
    return reverse_complement(my_seq)


def find_polyt_start(seq: str, window_size: int = 16, polya_fraction: float = 0.75) -> int:
    if len(seq) < window_size:
        return -1
    raw = np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)
    return int(find_polyt_start_bytes(raw, window_size, polya_fraction))


def find_polyt(seq: str, window_size: int = 16, polya_fraction: float = 0.75) -> Tuple[int, int]:
    polya_count = int(window_size * polya_fraction)
    n = len(seq)
    if n < window_size:
        return -1, -1
    i = 0
    a_count = seq[0:window_size].count("T")
    while i < n - window_size:
        if a_count >= polya_count:
            break
        first_t = seq[i] == "T"
        new_t = i + window_size < n and seq[i + window_size] == "T"
        if first_t and not new_t:
            a_count -= 1
        elif not first_t and new_t:
            a_count += 1
        i += 1

    if i >= n - window_size:
        return -1, -1
    polyt_start = i + max(0, seq[i:].find("TTTT"))

    while i < n - window_size:
        if a_count < polya_count:
            break
        first_t = seq[i] == "T"
        new_t = i + window_size < n and seq[i + window_size] == "T"
        if first_t and not new_t:
            a_count -= 1
        elif not first_t and new_t:
            a_count += 1
        i += 1

    if i >= n - window_size:
        polyt_end = n - 1
    else:
        last_t = seq[i : i + window_size].rfind("TTTT")
        polyt_end = i + window_size - 1 if last_t == -1 else i + last_t + 4

    return polyt_start, polyt_end


def align_pattern_ssw(
    sequence: str, start: int, end: int, pattern: str, min_score: int = 0
) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[int], Optional[int]]:
    seq = sequence[start:end]
    if not seq or not pattern:
        return None, None, None, None, None
    q_arr = np.frombuffer(pattern.encode("ascii"), dtype=np.uint8)
    r_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    score, ref_st, ref_en, read_st, read_en = sw_align_bytes(q_arr, r_arr, 1, 1, 1, 1)
    if score < min_score:
        return None, None, None, None, None
    return start + ref_st, start + ref_en, read_st, read_en, score


def find_candidate_with_max_score_ssw(
    barcode_matches: list,
    read_sequence: str,
    min_score: int = 10,
    score_diff: int = 0,
    sufficient_score: int = 0,
) -> Tuple[Optional[str], int, int, int]:
    if not read_sequence or not barcode_matches:
        return None, 0, 0, 0
    best_match = [0, 0, 0]
    best_barcode: Optional[str] = None
    second_best_score = 0

    r_arr = np.frombuffer(read_sequence.encode("ascii"), dtype=np.uint8)
    for barcode_match in barcode_matches:
        barcode = barcode_match[0]
        if not barcode:
            continue
        q_arr = np.frombuffer(barcode.encode("ascii"), dtype=np.uint8)
        score, ref_st, ref_en, read_st, read_en = sw_align_bytes(q_arr, r_arr, 1, 1, 1, 1)
        if score < min_score:
            continue

        if score > best_match[0]:
            best_barcode = barcode
            second_best_score = best_match[0]
            best_match[0] = score
            best_match[1] = ref_st - read_st
            best_match[2] = ref_en + (len(barcode) - read_en)
        elif score == best_match[0] and ref_st < best_match[1]:
            best_barcode = barcode
            second_best_score = best_match[0]
            best_match[1] = ref_st - read_st
            best_match[2] = ref_en + (len(barcode) - read_en)

        if score > sufficient_score > 0:
            break

    if best_match[0] - second_best_score < score_diff:
        return None, 0, 0, 0

    return best_barcode, best_match[0], best_match[1], best_match[2]


def find_candidate_with_max_score_ssw_var_len(
    barcode_matches: list, read_sequence: str, min_score: int = 14, score_diff: int = 1
) -> Tuple[Optional[str], int, int, int]:
    if not read_sequence or not barcode_matches:
        return None, 0, 0, 0
    best_match = [0, 0, 0, 0]
    second_best_match = [0, 0, 0, 0]
    best_barcode: Optional[str] = None

    r_arr = np.frombuffer(read_sequence.encode("ascii"), dtype=np.uint8)
    for barcode_match in barcode_matches:
        barcode = barcode_match[0]
        if not barcode:
            continue
        q_arr = np.frombuffer(barcode.encode("ascii"), dtype=np.uint8)
        score, ref_st, ref_en, read_st, read_en = sw_align_bytes(q_arr, r_arr, 1, 1, 1, 1)
        if score < min_score:
            continue

        ed = len(barcode) - score
        if score > best_match[0]:
            second_best_match = copy.copy(best_match)
            best_barcode = barcode
            best_match[0] = score
            best_match[1] = ref_st - read_st
            best_match[2] = ref_en + (len(barcode) - read_en)
            best_match[3] = ed
        elif score == best_match[0] and ed <= best_match[3]:
            second_best_match = copy.copy(best_match)
            best_barcode = barcode
            best_match[1] = ref_st - read_st
            best_match[2] = ref_en + (len(barcode) - read_en)
            best_match[3] = ed

    if best_barcode and best_match[0] < len(best_barcode) and best_match[0] - second_best_match[0] < score_diff:
        return None, 0, 0, 0

    return best_barcode, best_match[0], best_match[1], best_match[2]


def detect_exact_positions(
    sequence: str,
    start: int,
    end: int,
    kmer_size: int,
    pattern: str,
    pattern_occurrences: list,
    min_score: int = 0,
    start_delta: int = -1,
    end_delta: int = -1,
) -> Tuple[Optional[int], Optional[int]]:
    pattern_index = None
    for i, p in enumerate(pattern_occurrences):
        if p[0] == pattern:
            pattern_index = i
            break
    if pattern_index is None:
        return None, None

    start_pos = end_pos = pattern_start = pattern_end = None
    score = 0
    last_potential_pos = -2 * len(pattern)
    for match_position in pattern_occurrences[pattern_index][2]:
        if match_position - last_potential_pos < len(pattern):
            continue
        potential_start = max(start, start + match_position - len(pattern) + kmer_size)
        potential_end = min(end, start + match_position + len(pattern) + 1)
        al_r_st, al_r_en, al_q_st, al_q_en, al_sc = align_pattern_ssw(
            sequence, potential_start, potential_end, pattern, min_score
        )
        if al_sc is not None and al_sc > score:
            start_pos, end_pos, pattern_start, pattern_end, score = al_r_st, al_r_en, al_q_st, al_q_en, al_sc

    if start_pos is None or pattern_start is None or pattern_end is None or end_pos is None:
        return None, None
    if start_delta >= 0 and pattern_start > start_delta:
        return None, None
    if end_delta >= 0 and len(pattern) - pattern_end - 1 > end_delta:
        return None, None
    leftover_bases = len(pattern) - pattern_end - 1
    return max(0, start_pos - pattern_start), min(end_pos + leftover_bases, len(sequence) - 1)


def detect_first_exact_positions(
    sequence: str,
    start: int,
    end: int,
    kmer_size: int,
    pattern: str,
    pattern_occurrences: list,
    min_score: int = 0,
    start_delta: int = -1,
    end_delta: int = -1,
) -> Tuple[Optional[int], Optional[int]]:
    pattern_index = None
    for i, p in enumerate(pattern_occurrences):
        if p[0] == pattern:
            pattern_index = i
            break
    if pattern_index is None:
        return None, None

    start_pos = end_pos = pattern_start = pattern_end = None
    last_potential_pos = -2 * len(pattern)
    for match_position in pattern_occurrences[pattern_index][2]:
        if match_position - last_potential_pos < len(pattern):
            continue
        potential_start = max(start, start + match_position - len(pattern) + kmer_size)
        potential_end = min(end, start + match_position + len(pattern) + 1)
        al_r_st, al_r_en, al_q_st, al_q_en, al_sc = align_pattern_ssw(
            sequence, potential_start, potential_end, pattern, min_score
        )
        if al_sc is not None:
            start_pos, end_pos, pattern_start, pattern_end = al_r_st, al_r_en, al_q_st, al_q_en
            break

    if start_pos is None or pattern_start is None or pattern_end is None or end_pos is None:
        return None, None
    if start_delta >= 0 and pattern_start > start_delta:
        return None, None
    if end_delta >= 0 and len(pattern) - pattern_end - 1 > end_delta:
        return None, None
    leftover_bases = len(pattern) - pattern_end - 1
    return start_pos - pattern_start, end_pos + leftover_bases


def str_to_2bit(seq: str) -> int:
    kmer_idx = 0
    seq_len = len(seq)
    for i, ch in enumerate(seq):
        kmer_idx |= ((ord(ch) & 6) >> 1) << ((seq_len - i - 1) * 2)
    return kmer_idx


def bit_to_str(seq: int, seq_len: int) -> str:
    val = int(seq)
    return "".join(BIN2NUCL[(val >> ((seq_len - i - 1) * 2)) & 3] for i in range(seq_len))


def _convert_chunk_to_2bit(barcodes: List[str], seq_len: int) -> np.ndarray:
    n = len(barcodes)
    if n == 0:
        return np.array([], dtype=np.uint64)
    raw = "".join(barcodes).encode("ascii")
    byte_mat = np.ascontiguousarray(np.frombuffer(raw, dtype=np.uint8).reshape(n, seq_len))
    out = np.empty(n, dtype=np.uint64)
    encode_2bit_u8_to_u64(byte_mat, out)
    return out


def batch_str_to_2bit(barcodes: Iterable[str], seq_len: int) -> np.ndarray:
    bc_list = barcodes if isinstance(barcodes, list) else list(barcodes)
    if not bc_list:
        return np.array([], dtype=np.uint64)
    return _convert_chunk_to_2bit(bc_list, seq_len)


def batch_str_to_2bit_chunked(
    barcode_iterator: Iterable[str], seq_len: int, chunk_size: int = 10_000_000
) -> np.ndarray:
    chunks: List[np.ndarray] = []
    cur: List[str] = []
    for bc in barcode_iterator:
        cur.append(bc)
        if len(cur) >= chunk_size:
            chunks.append(_convert_chunk_to_2bit(cur, seq_len))
            cur.clear()
    if cur:
        chunks.append(_convert_chunk_to_2bit(cur, seq_len))
    if not chunks:
        return np.array([], dtype=np.uint64)
    return np.concatenate(chunks) if len(chunks) > 1 else chunks[0]


def edit_distance_fast(s1: str, s2: str) -> int:
    """Fast exact Cython Levenshtein distance replacing external editdistance package."""
    if s1 == s2:
        return 0
    b1 = np.frombuffer(s1.encode("ascii", errors="replace"), dtype=np.uint8)
    b2 = np.frombuffer(s2.encode("ascii", errors="replace"), dtype=np.uint8)
    return int(levenshtein_bytes(b1, b2))


def find_optimal_kmer_size(
    bc_len: int, bc_count: int, min_k: int = 6, max_k: int = 14, kmer_density: int = 50
) -> int:
    k = min_k
    while bc_count * (bc_len - k + 1) * (bc_len - k + 1) > (1 << (2 * k)) * kmer_density:
        k += 1
    return min(k, max_k)


def load_barcodes(inf: str | Path) -> List[str]:
    p = Path(inf)
    if not p.is_file():
        logger.critical("Barcode file '%s' does not exist.", p)
        sys.exit(int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND))
    opener = gzip.open if p.suffix in (".gz", ".gzip") else open
    with opener(p, "rt", encoding="utf-8") as fh:
        return [line.strip().split()[0] for line in fh if line.strip()]


def iter_barcodes(inf: str | Path) -> Iterator[str]:
    p = Path(inf)
    if not p.is_file():
        logger.critical("Barcode file '%s' does not exist.", p)
        sys.exit(int(IsoQuantExitCode.INPUT_FILE_NOT_FOUND))
    opener = gzip.open if p.suffix in (".gz", ".gzip") else open
    with opener(p, "rt", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if s:
                yield s.split()[0]


def _fallback_fastx_stream(path: Path) -> Iterator[Tuple[str, str, Optional[str]]]:
    opener = gzip.open if path.suffix in (".gz", ".gzip") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        first = ""
        for line in fh:
            if line.strip():
                first = line.rstrip("\r\n")
                break
        if not first:
            return
        if first.startswith("@"):
            # Streaming 4-line FASTQ fallback when synthetic files have mismatched quality lengths
            header: Optional[str] = first
            while header is not None:
                seq = next(fh, "").rstrip("\r\n")
                _plus = next(fh, "")
                qual = next(fh, "").rstrip("\r\n")
                rid = header[1:].split(None, 1)[0]
                yield rid, seq, qual
                header = None
                for nxt in fh:
                    if nxt.strip():
                        header = nxt.rstrip("\r\n")
                        break
        elif first.startswith(">"):
            cur_id = first[1:].split(None, 1)[0]
            chunks: List[str] = []
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                if s.startswith(">"):
                    yield cur_id, "".join(chunks), None
                    cur_id = s[1:].split(None, 1)[0]
                    chunks.clear()
                else:
                    chunks.append(s)
            yield cur_id, "".join(chunks), None


def iter_fastx_chunks(path: str | Path, chunk_size: int = 1000) -> Iterator[List[Tuple[str, str, Optional[str]]]]:
    """Stream (name, sequence, quality) batches via C/HTSlib pysam.FastxFile (no Bio.SeqRecord bloat)."""
    p = Path(path)
    yielded_ids = 0
    try:
        with pysam.FastxFile(str(p)) as fh:
            chunk: List[Tuple[str, str, Optional[str]]] = []
            for entry in fh:
                chunk.append((entry.name, entry.sequence or "", entry.quality))
                yielded_ids += 1
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
            if chunk:
                yield chunk
        return
    except ValueError:
        pass

    chunk = []
    for idx, item in enumerate(_fallback_fastx_stream(p)):
        if idx < yielded_ids:
            continue
        chunk.append(item)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
