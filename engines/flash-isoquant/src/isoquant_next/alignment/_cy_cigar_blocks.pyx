# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
"""Cython nogil CIGAR block extraction, aligned-pair error counting, and SoA locus buffer."""

from __future__ import annotations

import numpy as np
cimport numpy as cnp
from libc.stdint cimport int32_t, int64_t, uint8_t

cnp.import_array()


cpdef tuple extract_cigar_blocks_cy(
    const int32_t[:] ops,
    const int32_t[:] lens,
    int64_t ref_start,
):
    """
    Extract 1-based inclusive ref exons, 0-based inclusive read blocks, and cigar index blocks
    from CIGAR op/len arrays with a pure C nogil scan.
    Matches isoquant_lib.common.get_read_blocks identically.
    """
    cdef Py_ssize_t n = ops.shape[0]
    if n == 0:
        return [], [], []

    cdef cnp.ndarray[int64_t, ndim=1] ref_s = np.empty(n, dtype=np.int64)
    cdef cnp.ndarray[int64_t, ndim=1] ref_e = np.empty(n, dtype=np.int64)
    cdef cnp.ndarray[int64_t, ndim=1] rd_s = np.empty(n, dtype=np.int64)
    cdef cnp.ndarray[int64_t, ndim=1] rd_e = np.empty(n, dtype=np.int64)
    cdef cnp.ndarray[int32_t, ndim=1] cg_s = np.empty(n, dtype=np.int32)
    cdef cnp.ndarray[int32_t, ndim=1] cg_e = np.empty(n, dtype=np.int32)

    cdef int64_t[:] v_refs = ref_s
    cdef int64_t[:] v_refe = ref_e
    cdef int64_t[:] v_rds = rd_s
    cdef int64_t[:] v_rde = rd_e
    cdef int32_t[:] v_cgs = cg_s
    cdef int32_t[:] v_cge = cg_e

    cdef Py_ssize_t i, count = 0
    cdef int32_t op, ev_len
    cdef int64_t read_pos = 0
    cdef int64_t ref_pos = ref_start + 1
    cdef int64_t cur_ref_s = -1
    cdef int64_t cur_rd_s = -1
    cdef int32_t cur_cg_s = -1
    cdef bint has_block = 0
    cdef bint has_match = 0

    with nogil:
        for i in range(n):
            op = ops[i]
            ev_len = lens[i]
            if not has_block and (op == 0 or op == 1 or op == 2 or op == 7 or op == 8):
                has_block = 1
                cur_ref_s = ref_pos
                cur_rd_s = read_pos
                cur_cg_s = <int32_t>i
                if op == 1:
                    read_pos += ev_len
                elif op == 2:
                    ref_pos += ev_len
                else:
                    read_pos += ev_len
                    ref_pos += ev_len
                    has_match = 1
            elif op == 0 or op == 7 or op == 8:
                read_pos += ev_len
                ref_pos += ev_len
                has_match = 1
            elif op == 1:
                read_pos += ev_len
            elif op == 2:
                ref_pos += ev_len
            elif op == 3 or op == 4:
                if has_block:
                    if has_match:
                        v_refs[count] = cur_ref_s
                        v_refe[count] = ref_pos - 1
                        v_rds[count] = cur_rd_s
                        v_rde[count] = read_pos - 1
                        v_cgs[count] = cur_cg_s
                        v_cge[count] = <int32_t>(i - 1)
                        count += 1
                    has_block = 0
                    has_match = 0
                if op == 3:
                    ref_pos += ev_len
                else:
                    read_pos += ev_len

        if has_block and has_match:
            v_refs[count] = cur_ref_s
            v_refe[count] = ref_pos - 1
            v_rds[count] = cur_rd_s
            v_rde[count] = read_pos - 1
            v_cgs[count] = cur_cg_s
            v_cge[count] = <int32_t>(n - 1)
            count += 1

    cdef list ref_blocks = []
    cdef list read_blocks = []
    cdef list cigar_blocks = []
    for i in range(count):
        ref_blocks.append((int(v_refs[i]), int(v_refe[i])))
        read_blocks.append((int(v_rds[i]), int(v_rde[i])))
        cigar_blocks.append((int(v_cgs[i]), int(v_cge[i])))
    return ref_blocks, read_blocks, cigar_blocks


cpdef tuple count_errors_in_arrays_cy(
    const int32_t[:] q_pos,
    const int64_t[:] r_pos,
    int64_t ref_start,
    int64_t ref_end,
    const uint8_t[:] query_bytes,
    const uint8_t[:] ref_bytes,
    int64_t ref_offset,
):
    """
    Count (indel_count, mismatch_count) over pre-extracted aligned_pairs arrays
    matching AlignmentInfo.get_error_count 1-based coordinate rules in nogil C.
    Sentinels: q_pos == -1 means query_pos is None; r_pos == -1 means ref_pos is None.
    """
    cdef Py_ssize_t n = q_pos.shape[0]
    cdef Py_ssize_t i
    cdef int32_t qp
    cdef int64_t rp, true_rp
    cdef int64_t last_ref_pos = 0
    cdef int64_t indels = 0
    cdef int64_t mismatches = 0
    cdef bint check_seq = (query_bytes is not None and ref_bytes is not None)
    cdef Py_ssize_t qlen = query_bytes.shape[0] if query_bytes is not None else 0
    cdef Py_ssize_t rlen = ref_bytes.shape[0] if ref_bytes is not None else 0
    cdef int64_t chr_idx

    with nogil:
        for i in range(n):
            qp = q_pos[i]
            rp = r_pos[i]
            if rp < 0:
                if ref_start - 1 <= last_ref_pos <= ref_end:
                    indels += 1
            else:
                true_rp = rp + 1
                if ref_start <= true_rp <= ref_end:
                    if qp < 0:
                        indels += 1
                    elif check_seq and qp < qlen:
                        chr_idx = rp - ref_offset
                        if 0 <= chr_idx < rlen and query_bytes[qp] != ref_bytes[chr_idx]:
                            mismatches += 1
                last_ref_pos = true_rp
                if last_ref_pos > ref_end:
                    break
    return int(indels), int(mismatches)
