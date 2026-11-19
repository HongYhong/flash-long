# cython: language_level=3
"""Cython nogil kernels for local affine Smith-Waterman, Levenshtein UMI distance, 2-bit encoding, and polyT search."""

from libc.stdint cimport int8_t, int32_t, uint8_t, uint64_t
from libc.stdlib cimport free, malloc


cdef int _sw_align_c(
    const uint8_t* r,
    int m,
    const uint8_t* ref,
    int n,
    int match_score,
    int mismatch_penalty,
    int gap_open,
    int gap_ext,
    int* out_ref_start,
    int* out_ref_end,
    int* out_read_start,
    int* out_read_end,
) noexcept nogil:
    """Affine-gap local Smith-Waterman returning optimal score & 0-based inclusive endpoints."""
    if m <= 0 or n <= 0:
        out_ref_start[0] = 0
        out_ref_end[0] = -1
        out_read_start[0] = 0
        out_read_end[0] = -1
        return 0

    cdef int cols = n + 1
    cdef int* H_prev = <int*>malloc(cols * sizeof(int))
    cdef int* H_curr = <int*>malloc(cols * sizeof(int))
    cdef int* E_row = <int*>malloc(cols * sizeof(int))
    cdef uint8_t* trace = <uint8_t*>malloc((m + 1) * cols * sizeof(uint8_t))
    if H_prev == NULL or H_curr == NULL or E_row == NULL or trace == NULL:
        if H_prev != NULL: free(H_prev)
        if H_curr != NULL: free(H_curr)
        if E_row != NULL: free(E_row)
        if trace != NULL: free(trace)
        return -1

    # trace bits: low 2 bits = best state of H[i,j] (0=stop,1=diag,2=up/E,3=left/F)
    # bit 2 (4): E extended from E[i-1,j] (vs opened from H[i-1,j])
    # bit 3 (8): F extended from F[i,j-1] (vs opened from H[i,j-1])
    cdef int neg_inf = -1000000000
    cdef int i, j
    for j in range(cols):
        H_prev[j] = 0
        E_row[j] = neg_inf
        trace[j] = 0

    cdef int best_score = 0
    cdef int best_i = 0
    cdef int best_j = 0
    cdef int F_val, diag_val, e_open, e_ext, f_open, f_ext, h_val
    cdef uint8_t tb_code, st
    cdef int* tmp_ptr

    for i in range(1, m + 1):
        H_curr[0] = 0
        trace[i * cols] = 0
        F_val = neg_inf
        for j in range(1, cols):
            if r[i - 1] == ref[j - 1]:
                diag_val = H_prev[j - 1] + match_score
            else:
                diag_val = H_prev[j - 1] - mismatch_penalty

            e_open = H_prev[j] - gap_open
            e_ext = E_row[j] - gap_ext
            if e_open >= e_ext:
                E_row[j] = e_open
                tb_code = 0
            else:
                E_row[j] = e_ext
                tb_code = 4

            f_open = H_curr[j - 1] - gap_open
            f_ext = F_val - gap_ext
            if f_open >= f_ext:
                F_val = f_open
            else:
                F_val = f_ext
                tb_code |= 8

            h_val = 0
            st = 0
            if diag_val > h_val:
                h_val = diag_val
                st = 1
            if E_row[j] > h_val:
                h_val = E_row[j]
                st = 2
            if F_val > h_val:
                h_val = F_val
                st = 3

            H_curr[j] = h_val
            trace[i * cols + j] = tb_code | st
            if h_val > best_score:
                best_score = h_val
                best_i = i
                best_j = j

        tmp_ptr = H_prev
        H_prev = H_curr
        H_curr = tmp_ptr

    if best_score <= 0:
        free(H_prev)
        free(H_curr)
        free(E_row)
        free(trace)
        out_ref_start[0] = 0
        out_ref_end[0] = -1
        out_read_start[0] = 0
        out_read_end[0] = -1
        return 0

    out_read_end[0] = best_i - 1
    out_ref_end[0] = best_j - 1

    # Traceback exact state machine back to local alignment start (st == 0)
    cdef int cur_i = best_i
    cdef int cur_j = best_j
    cdef int mode = 1  # 1=H, 2=E(up), 3=F(left)
    cdef int start_i = cur_i
    cdef int start_j = cur_j
    cdef uint8_t code

    while cur_i > 0 and cur_j > 0:
        code = trace[cur_i * cols + cur_j]
        if mode == 1:
            st = code & 3
            if st == 0:
                break
            elif st == 1:
                start_i = cur_i - 1
                start_j = cur_j - 1
                cur_i -= 1
                cur_j -= 1
            elif st == 2:
                mode = 2
            else:
                mode = 3
        elif mode == 2:
            start_i = cur_i - 1
            start_j = cur_j - 1 if cur_j > 0 else 0
            if (code & 4) == 0:
                mode = 1
            cur_i -= 1
        else:
            start_i = cur_i - 1 if cur_i > 0 else 0
            start_j = cur_j - 1
            if (code & 8) == 0:
                mode = 1
            cur_j -= 1

    out_read_start[0] = start_i
    out_ref_start[0] = start_j

    free(H_prev)
    free(H_curr)
    free(E_row)
    free(trace)
    return best_score


def sw_align_bytes(
    const uint8_t[::1] read_bytes,
    const uint8_t[::1] ref_bytes,
    int match_score=1,
    int mismatch_penalty=1,
    int gap_open=1,
    int gap_extension=1,
):
    """Return (optimal_score, ref_start, ref_end, read_start, read_end) 0-based inclusive."""
    cdef int m = <int>read_bytes.shape[0]
    cdef int n = <int>ref_bytes.shape[0]
    if m == 0 or n == 0:
        return 0, 0, -1, 0, -1
    cdef int r_st = 0, r_en = -1, q_st = 0, q_en = -1, score = 0
    with nogil:
        score = _sw_align_c(
            &read_bytes[0], m, &ref_bytes[0], n,
            match_score, mismatch_penalty, gap_open, gap_extension,
            &r_st, &r_en, &q_st, &q_en,
        )
    if score < 0:
        raise MemoryError("Failed allocating SW DP buffers")
    return score, r_st, r_en, q_st, q_en


def levenshtein_bytes(const uint8_t[::1] a, const uint8_t[::1] b) -> int:
    """Exact Levenshtein edit distance in Cython nogil (1D rolling array)."""
    cdef int la = <int>a.shape[0]
    cdef int lb = <int>b.shape[0]
    if la == 0:
        return lb
    if lb == 0:
        return la
    if la > lb:
        return levenshtein_bytes(b, a)

    cdef int i, j, prev_diag, temp, cost
    cdef int[256] stack_row
    cdef int* row = stack_row
    cdef int heap_alloc = 0
    if la + 1 > 256:
        row = <int*>malloc((la + 1) * sizeof(int))
        if row == NULL:
            raise MemoryError()
        heap_alloc = 1

    with nogil:
        for i in range(la + 1):
            row[i] = i
        for j in range(1, lb + 1):
            prev_diag = row[0]
            row[0] = j
            for i in range(1, la + 1):
                temp = row[i]
                cost = 0 if a[i - 1] == b[j - 1] else 1
                if row[i] + 1 < row[i - 1] + 1:
                    row[i] = row[i] + 1
                else:
                    row[i] = row[i - 1] + 1
                if prev_diag + cost < row[i]:
                    row[i] = prev_diag + cost
                prev_diag = temp
    cdef int ans = row[la]
    if heap_alloc:
        free(row)
    return ans


def encode_2bit_u8_to_u64(const uint8_t[:, ::1] encoded_chars, uint64_t[::1] out_codes) -> None:
    """Encode N x L matrix of ASCII bases into uint64 2-bit packed integers (A=0,C=1,T=2,G=3)."""
    cdef Py_ssize_t n = encoded_chars.shape[0]
    cdef Py_ssize_t seq_len = encoded_chars.shape[1]
    cdef Py_ssize_t i, pos
    cdef uint64_t acc
    with nogil:
        for i in range(n):
            acc = 0
            for pos in range(seq_len):
                acc = (acc << 2) | ((<uint64_t>encoded_chars[i, pos] & <uint64_t>6) >> 1)
            out_codes[i] = acc


def find_polyt_start_bytes(
    const uint8_t[::1] seq,
    int window_size=16,
    double polya_fraction=0.75,
) -> int:
    """Exact streaming Cython equivalent of IsoQuant find_polyt_start."""
    cdef int n = <int>seq.shape[0]
    if n < window_size or window_size <= 0:
        return -1
    cdef int polya_count = <int>(window_size * polya_fraction)
    cdef int i = 0
    cdef int count_t = 0
    cdef int limit = n - window_size
    cdef int j
    cdef int ans = -1

    with nogil:
        for j in range(window_size):
            if seq[j] == 84:  # 'T'
                count_t += 1
        while i < limit:
            if count_t >= polya_count:
                break
            if seq[i] == 84 and (i + window_size >= n or seq[i + window_size] != 84):
                count_t -= 1
            elif seq[i] != 84 and (i + window_size < n and seq[i + window_size] == 84):
                count_t += 1
            i += 1
        if i < limit:
            ans = i
            for j in range(i, n - 3):
                if seq[j] == 84 and seq[j + 1] == 84 and seq[j + 2] == 84 and seq[j + 3] == 84:
                    ans = j
                    break
    return ans
