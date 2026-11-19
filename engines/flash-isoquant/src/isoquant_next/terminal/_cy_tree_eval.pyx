# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
"""L1-Cache friendly nogil float32 GBDT inference kernel (zero xgboost/pandas deps)."""

from libc.math cimport exp, isnan
from libc.stdint cimport int32_t, uint8_t
import numpy as np
cimport numpy as cnp


cpdef cnp.ndarray[cnp.float64_t, ndim=1] eval_gbtree_proba_cy(
    const float[:, :] X,
    const int32_t[:] tree_roots,
    const int32_t[:] left_child,
    const int32_t[:] right_child,
    const int32_t[:] feature_idx,
    const uint8_t[:] default_left,
    const float[:] split_val,
    double const_margin,
):
    """Evaluate packed XGBoost JSON gbtree ensemble over float32 feature matrix X in nogil C."""
    cdef Py_ssize_t n_rows = X.shape[0]
    cdef Py_ssize_t n_trees = tree_roots.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1] out = np.empty(n_rows, dtype=np.float64)
    cdef double[:] out_view = out

    cdef Py_ssize_t r, t_i
    cdef int32_t node, lc, fidx
    cdef float fval, thresh
    cdef double margin

    with nogil:
        for r in range(n_rows):
            margin = const_margin
            for t_i in range(n_trees):
                node = tree_roots[t_i]
                while True:
                    lc = left_child[node]
                    if lc == -1:
                        margin += <double>split_val[node]
                        break
                    fidx = feature_idx[node]
                    fval = X[r, fidx]
                    if isnan(fval):
                        if default_left[node] != 0:
                            node = lc
                        else:
                            node = right_child[node]
                    else:
                        thresh = split_val[node]
                        if fval < thresh:
                            node = lc
                        else:
                            node = right_child[node]
            if margin > 500.0:
                out_view[r] = 1.0
            elif margin < -500.0:
                out_view[r] = 0.0
            else:
                out_view[r] = 1.0 / (1.0 + exp(-margin))

    return out
