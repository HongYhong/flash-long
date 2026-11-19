# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
"""Cython nogil sparse (pos, count) histogram moments and padded buffer builder."""

from libc.math cimport log, sqrt, NAN
from libc.stdint cimport int32_t, int64_t
import numpy as np
cimport numpy as cnp


cpdef tuple sparse_hist_moments_and_pad_cy(
    const int64_t[:] u_pos,
    const int64_t[:] u_cnt,
    Py_ssize_t pad_len,
):
    """Compute exact (data_min, total_n, mode_off, var, skew, entropy, mean_height, padded) in one pass.

    Eliminates upstream ``[pos] * count`` list expansion & Python np.histogram loops.
    """
    cdef Py_ssize_t k = u_pos.shape[0]
    cdef int64_t data_min = u_pos[0]
    cdef int64_t data_max = u_pos[k - 1]
    cdef Py_ssize_t i
    cdef int64_t total_n = 0
    cdef int64_t best_cnt = -1
    cdef int64_t best_pos = data_min
    cdef double sum_x = 0.0
    cdef double sum_c_log_c = 0.0
    cdef double d, c_val

    with nogil:
        for i in range(k):
            c_val = <double>u_cnt[i]
            if u_cnt[i] <= 0:
                continue
            total_n += u_cnt[i]
            d = <double>(u_pos[i] - data_min)
            sum_x += c_val * d
            if u_cnt[i] > best_cnt:
                best_cnt = u_cnt[i]
                best_pos = u_pos[i]
            sum_c_log_c += c_val * log(c_val)

    cdef double fn = <double>total_n
    cdef double mu = sum_x / fn
    cdef double m2 = 0.0
    cdef double m3 = 0.0
    cdef double diff, diff2
    with nogil:
        for i in range(k):
            if u_cnt[i] <= 0:
                continue
            diff = <double>(u_pos[i] - data_min) - mu
            diff2 = diff * diff
            c_val = <double>u_cnt[i]
            m2 += c_val * diff2
            m3 += c_val * diff2 * diff
        m2 /= fn
        m3 /= fn

    cdef double skew_val = <double>NAN
    if total_n > 3 and m2 >= 1e-12:
        skew_val = m3 / (m2 * sqrt(m2))

    cdef double entropy_val = 0.0
    if total_n > 0:
        entropy_val = log(fn) - (sum_c_log_c / fn)
        if entropy_val < 0.0:
            entropy_val = 0.0

    cdef Py_ssize_t span = <Py_ssize_t>(data_max - data_min + 1)
    cdef Py_ssize_t padded_len = span + 2 * pad_len
    cdef cnp.ndarray[cnp.int32_t, ndim=1] padded = np.zeros(padded_len, dtype=np.int32)
    cdef int32_t[:] p_view = padded
    with nogil:
        for i in range(k):
            if u_cnt[i] > 0:
                p_view[pad_len + <Py_ssize_t>(u_pos[i] - data_min)] = <int32_t>u_cnt[i]

    cdef double mean_height = fn / <double>padded_len
    return (
        <int64_t>data_min,
        <int64_t>total_n,
        <int64_t>(best_pos - data_min),
        m2,
        skew_val,
        entropy_val,
        mean_height,
        padded,
    )
