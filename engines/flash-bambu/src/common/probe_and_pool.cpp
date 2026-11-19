#include <RcppArmadillo.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <thread>
#include <atomic>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include "S4Vectors_interface.h"
#include "IRanges_interface.h"
}

extern "C" SEXP bambu_next_probe_abi(SEXP ir_test_sexp, SEXP n_threads_sexp) {
    int n_threads = Rf_asInteger(n_threads_sexp);
    if (n_threads < 1) n_threads = 1;

    int ir_len = -1;
    if (!Rf_isNull(ir_test_sexp)) {
        ir_len = get_IRanges_length(ir_test_sexp);
    }

    // Smoke-test std::thread workers without touching any R API inside workers
    std::atomic<int> acc{0};
    int k = std::min(n_threads, 4);
    std::vector<std::thread> workers;
    workers.reserve(k);
    for (int i = 0; i < k; ++i) {
        workers.emplace_back([&acc, i]() {
            acc.fetch_add(i + 1, std::memory_order_relaxed);
        });
    }
    for (auto &t : workers) t.join();

    int omp_max = 1;
#ifdef _OPENMP
    omp_max = omp_get_max_threads();
#endif

    const char *hts_ver = hts_version();

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, Rf_mkChar("htslib_version"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("iranges_c_len"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("thread_sum"));
    SET_STRING_ELT(nms, 3, Rf_mkChar("omp_max_threads"));
    SET_VECTOR_ELT(out, 0, Rf_mkString(hts_ver ? hts_ver : "unknown"));
    SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(ir_len));
    SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(acc.load()));
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(omp_max));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(2);
    return out;
}
