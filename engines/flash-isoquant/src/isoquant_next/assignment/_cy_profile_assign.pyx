# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
"""Cython nogil interval-profile constructor & matcher over compact int8 arrays and mmap .iqidx views."""

from __future__ import annotations

import numpy as np
cimport numpy as cnp
from libc.stdint cimport int8_t, int32_t, int64_t

cnp.import_array()


cdef inline bint _overlaps(int64_t s1, int64_t e1, int64_t s2, int64_t e2) noexcept nogil:
    return not (e1 < s2 or s1 > e2)


cdef inline bint _contains(int64_t bs, int64_t be, int64_t ss, int64_t se) noexcept nogil:
    return be >= se and bs <= ss


cdef inline bint _overlaps_at_least(int64_t s1, int64_t e1, int64_t s2, int64_t e2, int64_t delta) noexcept nogil:
    cdef int64_t ovlp1 = e1 - s2
    cdef int64_t ovlp2 = e2 - s1
    cdef int64_t d = delta - 1
    if ovlp1 < 0 or ovlp2 < 0:
        return 0
    if e1 < e2:
        return ovlp1 >= d or s1 >= s2
    return ovlp2 >= d or s1 <= s2


cdef inline bint _overlaps_at_least_when_overlap(int64_t s1, int64_t e1, int64_t s2, int64_t e2, int64_t delta) noexcept nogil:
    if e1 < e2:
        return s1 >= s2 or (e1 - s2 + 1) >= delta
    return s1 <= s2 or (e2 - s1 + 1) >= delta


cdef inline bint _equal_ranges(int64_t s1, int64_t e1, int64_t s2, int64_t e2, int64_t delta) noexcept nogil:
    cdef int64_t d0 = s1 - s2 if s1 >= s2 else s2 - s1
    cdef int64_t d1 = e1 - e2 if e1 >= e2 else e2 - e1
    return d0 <= delta and d1 <= delta


cdef inline int64_t _abs_i64(int64_t x) noexcept nogil:
    return x if x >= 0 else -x


cpdef tuple overlapping_profile_cy(
    const int32_t[:] ref_s,
    const int32_t[:] ref_e,
    const int32_t[:] rd_s,
    const int32_t[:] rd_e,
    int64_t gene_s,
    int64_t gene_e,
    int64_t map_s,
    int64_t map_e,
    int32_t comp_delta,
    int32_t abs_kind,
    int32_t abs_delta,
    int64_t polya_pos,
    int64_t polyt_pos,
    int32_t tail_delta,
):
    """
    Construct compact int8 gene/read profiles for possibly overlapping known features (introns/exons)
    without allocating Python list[int] or dict objects inside the scan loop.
    abs_kind: 0 -> contains, 1 -> overlaps_at_least(abs_delta)
    """
    cdef Py_ssize_t n_ref = ref_s.shape[0]
    cdef Py_ssize_t n_rd = rd_s.shape[0]

    cdef cnp.ndarray[int8_t, ndim=1] g_prof_arr = np.zeros(n_ref, dtype=np.int8)
    cdef cnp.ndarray[int8_t, ndim=1] r_prof_arr = np.zeros(n_rd, dtype=np.int8)
    # Track first matched gene index & best delta per read feature for duplicate tie-break in C
    cdef cnp.ndarray[int32_t, ndim=1] best_g_idx = np.full(n_rd, -1, dtype=np.int32)
    cdef cnp.ndarray[int32_t, ndim=1] first_g_idx = np.full(n_rd, -1, dtype=np.int32)
    cdef cnp.ndarray[int64_t, ndim=1] best_delta = np.full(n_rd, (<int64_t>1) << 60, dtype=np.int64)
    cdef cnp.ndarray[int32_t, ndim=1] match_count = np.zeros(n_rd, dtype=np.int32)

    cdef int8_t[:] gp = g_prof_arr
    cdef int8_t[:] rp = r_prof_arr
    cdef int32_t[:] bg = best_g_idx
    cdef int32_t[:] fg = first_g_idx
    cdef int64_t[:] bd = best_delta
    cdef int32_t[:] mc = match_count

    cdef Py_ssize_t i, gene_pos = 0, read_pos = 0
    cdef int64_t rs, re, is_, ie, d
    cdef bint cond_abs

    with nogil:
        for i in range(n_ref):
            if abs_kind == 0:
                cond_abs = _contains(map_s, map_e, ref_s[i], ref_e[i])
            else:
                cond_abs = _overlaps_at_least(map_s, map_e, ref_s[i], ref_e[i], abs_delta)
            if cond_abs:
                gp[i] = -1

        for i in range(n_rd):
            if abs_kind == 0:
                cond_abs = _contains(gene_s, gene_e, rd_s[i], rd_e[i])
            else:
                cond_abs = _overlaps_at_least(gene_s, gene_e, rd_s[i], rd_e[i], abs_delta)
            if cond_abs:
                rp[i] = -1

        while gene_pos < n_ref and read_pos < n_rd:
            rs = rd_s[read_pos]
            re = rd_e[read_pos]
            is_ = ref_s[gene_pos]
            ie = ref_e[gene_pos]
            if re < is_:
                if rp[read_pos] == 0 and gene_pos > 0:
                    rp[read_pos] = -1
                read_pos += 1
            elif ie < rs:
                if read_pos > 0:
                    gp[gene_pos] = -1
                gene_pos += 1
            elif _equal_ranges(rs, re, is_, ie, comp_delta):
                gp[gene_pos] = 1
                rp[read_pos] = 1
                if fg[read_pos] == -1:
                    fg[read_pos] = <int32_t>gene_pos
                d = _abs_i64(rs - is_) + _abs_i64(re - ie)
                mc[read_pos] += 1
                if d < bd[read_pos]:
                    bd[read_pos] = d
                    bg[read_pos] = <int32_t>gene_pos
                gene_pos += 1
            elif _overlaps(rs, re, is_, ie):
                if abs_kind == 0:
                    cond_abs = _contains(map_s, map_e, is_, ie)
                else:
                    cond_abs = _overlaps_at_least(map_s, map_e, is_, ie, abs_delta)
                if cond_abs:
                    gp[gene_pos] = -1
                gene_pos += 1

        # Second scan only if any read feature matched > 1 reference feature (eliminate suboptimal matches)
        for i in range(n_rd):
            if mc[i] > 1:
                rs = rd_s[i]
                re = rd_e[i]
                d = bd[i]
                for gene_pos in range(fg[i], n_ref):
                    is_ = ref_s[gene_pos]
                    ie = ref_e[gene_pos]
                    if ie < rs and not _equal_ranges(rs, re, is_, ie, comp_delta):
                        continue
                    if is_ > re + comp_delta:
                        break
                    if _equal_ranges(rs, re, is_, ie, comp_delta):
                        if _abs_i64(rs - is_) + _abs_i64(re - ie) > d:
                            gp[gene_pos] = -1

        if polya_pos != -1:
            for i in range(n_ref):
                if ref_s[i] > polya_pos + tail_delta:
                    gp[i] = -2
        if polyt_pos != -1:
            for i in range(n_ref):
                if ref_e[i] < polyt_pos - tail_delta:
                    gp[i] = -2

    cdef Py_ssize_t start_pos = 0
    cdef Py_ssize_t end_pos = n_ref - 1
    with nogil:
        while start_pos < n_ref and gp[start_pos] == 0:
            start_pos += 1
        while end_pos >= 0 and gp[end_pos] == 0:
            end_pos -= 1

    return g_prof_arr, r_prof_arr, first_g_idx, (<Py_ssize_t>start_pos, <Py_ssize_t>(end_pos + 1))


cdef inline Py_ssize_t _interval_bin_search(const int32_t[:] s_arr, const int32_t[:] e_arr, Py_ssize_t n, int64_t pos) noexcept nogil:
    if n == 0 or pos < s_arr[0] or pos > e_arr[n - 1]:
        return -1
    if pos >= s_arr[n - 1]:
        return n - 1
    cdef Py_ssize_t lo = 0
    cdef Py_ssize_t hi = n - 1
    cdef Py_ssize_t mid
    while lo < hi:
        mid = (lo + hi + 1) >> 1
        if s_arr[mid] <= pos:
            lo = mid
        else:
            hi = mid - 1
    return lo


cdef inline Py_ssize_t _interval_bin_search_rev(const int32_t[:] s_arr, const int32_t[:] e_arr, Py_ssize_t n, int64_t pos) noexcept nogil:
    if n == 0 or pos < s_arr[0] or pos > e_arr[n - 1]:
        return -1
    if pos <= e_arr[0]:
        return 0
    cdef Py_ssize_t lo = 0
    cdef Py_ssize_t hi = n - 1
    cdef Py_ssize_t mid
    while lo < hi:
        mid = (lo + hi) >> 1
        if e_arr[mid] >= pos:
            hi = mid
        else:
            lo = mid + 1
    return lo


cpdef tuple nonoverlapping_profile_cy(
    const int32_t[:] ref_s,
    const int32_t[:] ref_e,
    const int32_t[:] rd_s,
    const int32_t[:] rd_e,
    int32_t comp_kind,
    int32_t comp_delta,
    int32_t delta,
    int64_t polya_pos,
    int64_t polyt_pos,
):
    """
    Construct compact int8 profiles over disjoint reference split exons in nogil C.
    comp_kind: 0 -> overlaps, 1 -> overlaps_at_least_when_overlap(comp_delta)
    """
    cdef Py_ssize_t n_ref = ref_s.shape[0]
    cdef Py_ssize_t n_rd = rd_s.shape[0]
    cdef cnp.ndarray[int8_t, ndim=1] g_prof_arr = np.zeros(n_ref, dtype=np.int8)
    cdef cnp.ndarray[int8_t, ndim=1] r_prof_arr = np.zeros(n_rd, dtype=np.int8)
    cdef int8_t[:] ep = g_prof_arr
    cdef int8_t[:] rp = r_prof_arr

    cdef Py_ssize_t gene_pos = 0, read_pos = 0, i, pidx
    cdef int64_t rs, re, gs, ge
    cdef bint ok

    with nogil:
        while gene_pos < n_ref and read_pos < n_rd:
            rs = rd_s[read_pos]
            re = rd_e[read_pos]
            gs = ref_s[gene_pos]
            ge = ref_e[gene_pos]
            if re < gs:
                if gene_pos > 0 and rp[read_pos] == 0:
                    rp[read_pos] = -1
                read_pos += 1
            elif ge < rs:
                if read_pos > 0 and ep[gene_pos] == 0:
                    ep[gene_pos] = -1
                gene_pos += 1
            else:
                if comp_kind == 0:
                    ok = 1
                else:
                    ok = _overlaps_at_least_when_overlap(rs, re, gs, ge, comp_delta)
                if ok:
                    ep[gene_pos] = 1
                    rp[read_pos] = 1
                if re < ge:
                    read_pos += 1
                else:
                    gene_pos += 1

        if polya_pos != -1:
            pidx = _interval_bin_search(ref_s, ref_e, n_ref, polya_pos + delta)
            if pidx != -1:
                for i in range(pidx + 1, n_ref):
                    ep[i] = -2
        if polyt_pos != -1:
            pidx = _interval_bin_search_rev(ref_s, ref_e, n_ref, polyt_pos - delta)
            if pidx != -1:
                for i in range(0, pidx):
                    ep[i] = -2

    cdef Py_ssize_t start_pos = 0
    cdef Py_ssize_t end_pos = n_ref - 1
    with nogil:
        while start_pos < n_ref and ep[start_pos] == 0:
            start_pos += 1
        while end_pos >= 0 and ep[end_pos] == 0:
            end_pos -= 1

    return g_prof_arr, r_prof_arr, (<Py_ssize_t>start_pos, <Py_ssize_t>(end_pos + 1))


cpdef Py_ssize_t diff_present_features_cy(
    const int8_t[:] prof1,
    const int8_t[:] prof2,
    Py_ssize_t s,
    Py_ssize_t e,
    Py_ssize_t diff_limit,
):
    cdef Py_ssize_t d = 0
    cdef Py_ssize_t i
    cdef int8_t p1, p2
    cdef Py_ssize_t limit = prof1.shape[0] + 1 if diff_limit == -1 else diff_limit
    with nogil:
        for i in range(s, e):
            p1 = prof1[i]
            p2 = prof2[i]
            if p1 == 0 or p2 == 0:
                continue
            if p1 != p2:
                d += 1
                if d > limit:
                    break
    return d


cpdef bint has_overlapping_features_cy(
    const int8_t[:] prof1,
    const int8_t[:] prof2,
    Py_ssize_t s,
    Py_ssize_t e,
):
    cdef Py_ssize_t i
    cdef bint found = 0
    with nogil:
        for i in range(s, e):
            if prof1[i] == 1 and prof2[i] == 1:
                found = 1
                break
    return found


cpdef int64_t count_iqidx_exon_overlaps_cy(
    const int32_t[:] iq_exon_s,
    const int32_t[:] iq_exon_e,
    Py_ssize_t tx_exon_start,
    Py_ssize_t tx_exon_end,
    const int32_t[:] rd_s,
    const int32_t[:] rd_e,
):
    """Directly scan mmap .iqidx exon slice against sorted read exons in nogil C without copying."""
    cdef Py_ssize_t i = tx_exon_start
    cdef Py_ssize_t j = 0
    cdef Py_ssize_t n_rd = rd_s.shape[0]
    cdef int64_t overlap_bp = 0
    cdef int32_t s1, e1, s2, e2
    with nogil:
        while i < tx_exon_end and j < n_rd:
            s1 = iq_exon_s[i]
            e1 = iq_exon_e[i]
            s2 = rd_s[j]
            e2 = rd_e[j]
            if not (e1 < s2 or s1 > e2):
                overlap_bp += (e1 if e1 < e2 else e2) - (s1 if s1 > s2 else s2) + 1
            if e1 < e2:
                i += 1
            else:
                j += 1
    return overlap_bp
