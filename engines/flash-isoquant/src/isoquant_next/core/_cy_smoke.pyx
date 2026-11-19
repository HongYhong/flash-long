# cython: language_level=3
"""Smoke Cython extension used to validate in-place build toolchain in sqanti3."""

from libc.stdint cimport int64_t


def fast_interval_sum(const int[::1] starts, const int[::1] ends) -> int64_t:
    """Compute total span across aligned coordinate pairs without GIL."""
    cdef Py_ssize_t n = starts.shape[0]
    cdef Py_ssize_t i
    cdef int64_t total = 0
    if n != ends.shape[0]:
        raise ValueError("starts and ends must have equal length")
    with nogil:
        for i in range(n):
            if ends[i] > starts[i]:
                total += ends[i] - starts[i]
    return total
