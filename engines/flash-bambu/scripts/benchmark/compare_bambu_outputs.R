#!/usr/bin/env Rscript
# Strict L1 (structure identical) and L2 (quantitative numeric bounds) comparator
suppressPackageStartupMessages({
    library(methods)
    library(SummarizedExperiment)
    library(GenomicRanges)
    library(S4Vectors)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: compare_bambu_outputs.R <mode:s1|s2|full> <ref_path> <cand_path> [out_json]")
}
cmp_mode <- args[1]
ref_path <- args[2]
cand_path <- args[3]
out_json <- if (length(args) >= 4) args[4] else ""

read_obj <- function(p) {
    if (grepl("\\.rds$", p, ignore.case = TRUE)) return(readRDS(p))
    if (grepl("\\.qs$", p, ignore.case = TRUE)) return(qs::qread(p))
    stop("Unsupported file format: ", p)
}

compare_s1 <- function(ref, cand) {
    if (identical(ref, cand)) {
        return(list(pass = TRUE, l1_identical = TRUE, details = "identical"))
    }
    rr_eq <- identical(rowRanges(ref), rowRanges(cand))
    rd_ref <- rowData(ref); rd_cand <- rowData(cand)
    cols_eq <- identical(colnames(rd_ref), colnames(rd_cand))
    diff_cols <- character()
    if (cols_eq) {
        for (nm in colnames(rd_ref)) {
            if (!identical(rd_ref[[nm]], rd_cand[[nm]])) {
                ae <- all.equal(rd_ref[[nm]], rd_cand[[nm]], tolerance = 1e-12)
                if (!isTRUE(ae)) diff_cols <- c(diff_cols, paste0(nm, ":", paste(ae, collapse=";")))
            }
        }
    }
    pass <- isTRUE(rr_eq) && isTRUE(cols_eq) && length(diff_cols) == 0 && identical(assays(ref), assays(cand))
    list(pass = pass, l1_identical = identical(ref, cand), rr_identical = rr_eq,
         cols_identical = cols_eq, diff_cols = diff_cols)
}

compare_full_se <- function(ref, cand) {
    # L1 checks on annotations & rowData metadata
    rr_identical <- identical(rowRanges(ref), rowRanges(cand))
    rd_ref <- rowData(ref); rd_cand <- rowData(cand)
    tx_id_eq <- identical(rownames(ref), rownames(cand)) &&
        identical(as.character(rd_ref$TXNAME), as.character(rd_cand$TXNAME)) &&
        identical(as.character(rd_ref$GENEID), as.character(rd_cand$GENEID))
    col_eq <- identical(colnames(ref), colnames(cand))

    # L2 checks on all assays
    assay_names <- names(assays(ref))
    max_abs_err <- list()
    pearson_r <- list()
    l2_pass <- identical(assay_names, names(assays(cand)))
    for (an in assay_names) {
        m1 <- as.matrix(assays(ref)[[an]])
        m2 <- as.matrix(assays(cand)[[an]])
        if (!identical(dim(m1), dim(m2))) {
            l2_pass <- FALSE
            max_abs_err[[an]] <- Inf
            pearson_r[[an]] <- NA_real_
        } else {
            diff_v <- abs(m1 - m2)
            mx <- max(diff_v)
            max_abs_err[[an]] <- mx
            v1 <- as.numeric(m1); v2 <- as.numeric(m2)
            r_val <- if (sd(v1) == 0 && sd(v2) == 0) 1.0 else suppressWarnings(cor(v1, v2))
            pearson_r[[an]] <- r_val
            if (!is.finite(mx) || mx > 1e-5 || !is.finite(r_val) || r_val < 0.999999) {
                l2_pass <- FALSE
            }
        }
    }
    # Check incompatibleCounts if present
    inc_ref <- metadata(ref)$incompatibleCounts
    inc_cand <- metadata(cand)$incompatibleCounts
    inc_eq <- identical(inc_ref, inc_cand)
    if (!inc_eq && !is.null(inc_ref) && !is.null(inc_cand)) {
        inc_eq <- isTRUE(all.equal(inc_ref, inc_cand, tolerance = 1e-6))
    }

    list(
        pass = isTRUE(rr_identical && tx_id_eq && col_eq && l2_pass && inc_eq),
        strict_identical = identical(ref, cand),
        rowRanges_identical = rr_identical,
        tx_gene_identical = tx_id_eq,
        col_identical = col_eq,
        incompatible_equal = inc_eq,
        max_abs_err = max_abs_err,
        pearson_r = pearson_r
    )
}

ref_obj <- read_obj(ref_path)
cand_obj <- read_obj(cand_path)

res <- if (cmp_mode == "s1") {
    compare_s1(ref_obj, cand_obj)
} else if (cmp_mode == "full") {
    compare_full_se(ref_obj, cand_obj)
} else if (cmp_mode == "exact") {
    ae <- all.equal(ref_obj, cand_obj, tolerance = 1e-12)
    list(pass = identical(ref_obj, cand_obj) || isTRUE(ae),
         strict_identical = identical(ref_obj, cand_obj),
         all_equal = if (isTRUE(ae)) "TRUE" else paste(head(ae, 10), collapse = " | "))
} else {
    stop("Unknown cmp_mode: ", cmp_mode)
}

json_txt <- jsonlite::toJSON(res, auto_unbox = TRUE, pretty = TRUE)
cat(json_txt, "\n")
if (nzchar(out_json)) writeLines(json_txt, out_json)
if (!isTRUE(res$pass)) quit(status = 1L)
