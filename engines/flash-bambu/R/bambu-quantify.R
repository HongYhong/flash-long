#' Perform quantification
#' @inheritParams bambu
#' @import data.table
#' @noRd
bambu.quantify <- function(readClass, annotations, emParameters, 
                           trackReads = FALSE, returnDistTable = FALSE,
                           verbose = FALSE, isoreParameters = setIsoreParameters(NULL),
                           ncore = 1L) {
    start.ptm <- proc.time()
    if (is.character(readClass)) readClass <- readRDS(file = readClass)
    colNameRC <- colnames(readClass)
    colDataRC <- colData(readClass)
    rr <- SummarizedExperiment::rowRanges(readClass)
    n_q <- length(rr)

    if (n_q > 0L) {
        ex_cnts <- as.integer(IRanges::elementNROWS(rr))
        unlisted <- unlist(rr, use.names = FALSE)
        ex_mc <- S4Vectors::mcols(unlisted, use.names = FALSE)
        st_chr <- as.character(GenomicRanges::strand(unlisted))
        q_seq_names <- as.character(Seqinfo::seqlevels(rr))
        q_ex_seq_1b <- as.integer(GenomicRanges::seqnames(unlisted))
        q_ex_start <- as.integer(IRanges::start(unlisted))
        q_ex_end <- as.integer(IRanges::end(unlisted))
        q_ex_strand <- ifelse(st_chr == "+", 1L, ifelse(st_chr == "-", -1L, 0L))
        q_ex_rank <- if (!is.null(ex_mc$exon_rank)) as.integer(ex_mc$exon_rank) else seq_along(q_ex_start)
        q_ex_end_rank <- if (!is.null(ex_mc$exon_endRank)) as.integer(ex_mc$exon_endRank) else q_ex_rank
        q_read_count <- as.numeric(SummarizedExperiment::assays(readClass)$counts[, 1L])
        q_rc_names <- as.character(rownames(readClass))
    } else {
        q_seq_names <- character(0)
        ex_cnts <- integer(0)
        q_ex_seq_1b <- integer(0)
        q_ex_start <- integer(0)
        q_ex_end <- integer(0)
        q_ex_strand <- integer(0)
        q_ex_rank <- integer(0)
        q_ex_end_rank <- integer(0)
        q_read_count <- numeric(0)
        q_rc_names <- character(0)
    }

    csr_ptr <- bambu_next_get_csr(annotations, n_threads = ncore)
    iso_p <- list(
        min_exonDistance = as.integer(isoreParameters[["min.exonDistance"]]),
        min_primarySecondaryDist = as.integer(isoreParameters[["min.primarySecondaryDist"]]),
        min_primarySecondaryDistStartEnd2 = as.integer(isoreParameters[["min.primarySecondaryDistStartEnd2"]])
    )
    em_p <- list(
        degradationBias = isTRUE(emParameters[["degradationBias"]]),
        maxiter = as.integer(emParameters[["maxiter"]]),
        conv = as.numeric(emParameters[["conv"]]),
        minvalue = as.numeric(emParameters[["minvalue"]]),
        sig_digit = as.integer(emParameters[["sig.digit"]])
    )

    need_dt <- isTRUE(returnDistTable) || isTRUE(trackReads)
    res <- .Call(
        "bambu_next_quantify_cpp",
        csr_ptr,
        q_seq_names,
        ex_cnts,
        q_ex_seq_1b,
        q_ex_start,
        q_ex_end,
        q_ex_strand,
        q_ex_rank,
        q_ex_end_rank,
        q_read_count,
        q_rc_names,
        iso_p,
        em_p,
        need_dt,
        as.integer(ncore),
        PACKAGE = "bambu"
    )

    incomp_val <- res$incomp_counts
    if (length(incomp_val) > 0L && all(!is.na(incomp_val)) &&
        all(floor(incomp_val) == incomp_val) && max(incomp_val) <= 2147483647) {
        incomp_val <- as.integer(incomp_val)
    }
    incomp_dt <- data.table::data.table(GENEID = res$incomp_genes, counts = incomp_val)
    data.table::setnames(incomp_dt, "counts", colNameRC)

    seOutput <- SummarizedExperiment(
        assays = SimpleList(
            counts = matrix(res$counts, ncol = 1L, dimnames = list(NULL, colNameRC)),
            CPM = matrix(res$CPM, ncol = 1L, dimnames = list(NULL, colNameRC)),
            fullLengthCounts = matrix(res$fullLengthCounts, ncol = 1L, dimnames = list(NULL, colNameRC)),
            uniqueCounts = matrix(res$uniqueCounts, ncol = 1L, dimnames = list(NULL, colNameRC))
        ),
        colData = colDataRC
    )
    metadata(seOutput)$incompatibleCounts <- incomp_dt

    if (need_dt) {
        dt_raw <- res$dist_table
        rc_val <- dt_raw$readCount
        if (length(rc_val) > 0L && all(!is.na(rc_val)) &&
            all(floor(rc_val) == rc_val) && max(rc_val) <= 2147483647) {
            rc_val <- as.integer(rc_val)
        }
        dist_df <- S4Vectors::DataFrame(
            annotationTxId = dt_raw$annotationTxId,
            txid = dt_raw$txid,
            readClassId = dt_raw$readClassId,
            readCount = rc_val,
            compatible = dt_raw$compatible,
            equal = dt_raw$equal,
            dist = dt_raw$dist,
            GENEID = dt_raw$GENEID
        )
        if (isTRUE(returnDistTable)) {
            dist_df$eqClassById <- as.list(S4Vectors::mcols(annotations)$eqClassById)[dt_raw$txid]
            dist_df <- unique(dist_df)
            metadata(seOutput)$distTable <- dist_df
        } else {
            dist_df <- unique(dist_df)
        }
        if (isTRUE(trackReads)) {
            metadata(seOutput)$readToTranscriptMap <- generateReadToTranscriptMap(
                readClass, dist_df, annotations
            )
        }
    }
    end.ptm <- proc.time()
    if (verbose) message("Finished isoform quantification for ", colNameRC, " in ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    return(seOutput)
}



#' Process data.table object
#' @param readClassDt A data.table object
#' @inheritParams bambu
#' @noRd
bambu.quantDT <- function(readClassDt = readClassDt, 
                          emParameters = list(degradationBias = TRUE, maxiter = 10000, conv = 10^(-2),
                                              minvalue = 10^(-8)), ncore = 1, verbose = FALSE) {
    rcPreOut <- addAval(readClassDt, emParameters, verbose)
    readClassDt <- rcPreOut[[1]]
    outIni <- initialiseOutput(readClassDt)
    readClassDt <- filterTxRc(readClassDt) 
    readClassDt <- assignGroups(readClassDt)
    inputRcDt <- getInputList(readClassDt)
    readClassDt <- split(readClassDt, by = "gene_grp_id")
    start.ptm <- proc.time()
    outEst <- abundance_quantification(inputRcDt, readClassDt,
                                       maxiter = emParameters[["maxiter"]],
                                       conv = emParameters[["conv"]], minvalue = emParameters[["minvalue"]])
    end.ptm <- proc.time()
    if (verbose) message("Finished EM estimation in ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    outEst <- modifyQuantOut(outEst,outIni)
    theta_est <- rbind(rcPreOut[[2]],outEst)
    theta_est <- removeDuplicates(theta_est)
    return(theta_est)
}


