# Node 1 (S0) adapters keeping legacy function names for seamless compatibility

#' @noRd
prepareAnnotationsFromGTF <- function(file) {
    if (missing(file)) stop("A GTF file is required.")
    n_thr <- suppressWarnings(as.integer(Sys.getenv("BAMBU_NEXT_PREP_THREADS", unset = "4")))
    if (is.na(n_thr) || n_thr < 1L) n_thr <- 1L
    bambu_next_prepare_annotations_fast(file, n_threads = n_thr)
}

#' @noRd
getMinimumEqClassByTx <- function(exonsByTranscripts) {
    n_tx <- length(exonsByTranscripts)
    if (n_tx == 0L) {
        return(S4Vectors::DataFrame(queryTxId = character(0), eqClassById = IRanges::IntegerList()))
    }
    if (!identical(S4Vectors::mcols(exonsByTranscripts)$txid, seq_len(n_tx))) {
        return(S4Vectors::DataFrame(queryTxId = names(exonsByTranscripts), eqClassById = rep(NA, n_tx)))
    }
    n_thr <- suppressWarnings(as.integer(Sys.getenv("BAMBU_NEXT_PREP_THREADS", unset = "4")))
    if (is.na(n_thr) || n_thr < 1L) n_thr <- 1L
    # Force compute via C++ FlatCSR (ignoring stale eqClassById on exonsByTranscripts)
    mc_old <- S4Vectors::mcols(exonsByTranscripts)
    S4Vectors::mcols(exonsByTranscripts)$eqClassById <- NULL
    attr(exonsByTranscripts, "bambu_next_csr") <- NULL
    updated <- bambu_next_attach_csr(exonsByTranscripts, n_threads = n_thr, compute_eq = TRUE)
    S4Vectors::DataFrame(
        queryTxId = names(updated),
        eqClassById = S4Vectors::mcols(updated)$eqClassById
    )
}
