#' Prepare reference annotations (Node 1: C++ streaming GTF & Flat-CSR minEqClass)
#' @param x A path to gtf file or a \code{TxDb} object.
#' @return A \code{GRangesList} object
#' @export
prepareAnnotations <- function(x) {
    bambu_next_apply_thread_guardrails()
    n_thr <- suppressWarnings(as.integer(Sys.getenv("BAMBU_NEXT_PREP_THREADS", unset = "4")))
    if (is.na(n_thr) || n_thr < 1L) n_thr <- 1L
    if (!is(x, "TxDb")) {
        if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
            stop("Input annotation file not readable. Requires .gtf/.gff format or TxDb object")
        }
    }
    bambu_next_prepare_annotations_fast(x, n_threads = n_thr)
}
