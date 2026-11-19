# Node 0: Runtime guardrails, native ABI wrapper, and isolated BiocFileCache setup

#' @noRd
em_theta <- function(X, Y, maxiter, minvalue, conv) {
    .Call("bambu_next_em_theta", X, Y, as.integer(maxiter), as.numeric(minvalue), as.numeric(conv), PACKAGE = "bambu")
}

#' @noRd
emWithL1 <- function(A, A_full, A_unique, Y, K, maxiter, minvalue, conv) {
    .Call("bambu_next_emWithL1", A, A_full, A_unique, Y, K, as.integer(maxiter), as.numeric(minvalue), as.numeric(conv), PACKAGE = "bambu")
}

#' @noRd
bambu_next_apply_thread_guardrails <- function() {
    # Prevent nested multi-threading contention in data.table / xgboost / BLAS
    suppressWarnings({
        data.table::setDTthreads(1L)
        xgboost::xgb.set.config(nthread = 1L)
    })
    invisible(TRUE)
}

#' @noRd
bambu_next_init_isolated_bfc <- function(cache_dir = tempdir()) {
    # Initialize an isolated BiocFileCache in main process before any parallel fork
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    bfc_dir <- file.path(cache_dir, "bambu_next_bfc")
    if (!dir.exists(bfc_dir)) dir.create(bfc_dir, recursive = TRUE, showWarnings = FALSE)
    suppressMessages(BiocFileCache::BiocFileCache(cache = bfc_dir, ask = FALSE))
    bfc_dir
}

#' @noRd
bambu_next_info <- function(n_threads = 2L) {
    ir <- IRanges::IRanges(start = c(1L, 10L), end = c(5L, 20L))
    abi <- .Call("bambu_next_probe_abi", ir, as.integer(n_threads), PACKAGE = "bambu")
    list(
        package_version = as.character(utils::packageVersion("bambu")),
        implementation = "bambu_next",
        abi = abi
    )
}
