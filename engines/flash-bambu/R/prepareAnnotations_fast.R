# Node 1 (S0): Streaming C++ GTF -> GRangesList + FlatCSR + precomputed minEqClass

#' @noRd
.bambu_next_extract_digest_args <- function(annotations) {
    n_tx <- length(annotations)
    mc <- S4Vectors::mcols(annotations, use.names = FALSE)
    ex_starts <- annotations@unlistData@ranges@start
    ex_ends <- ex_starts + annotations@unlistData@ranges@width - 1L
    eq_lens <- NULL
    eq_vals <- NULL
    if (!is.null(mc$eqClassById)) {
        eq_list <- mc$eqClassById
        eq_lens <- as.integer(IRanges::elementNROWS(eq_list))
        eq_vals <- as.integer(unlist(eq_list, use.names = FALSE))
    }
    first_tx <- if (!is.null(mc$TXNAME)) as.character(mc$TXNAME[1L]) else ""
    last_tx <- if (!is.null(mc$TXNAME)) as.character(mc$TXNAME[n_tx]) else ""
    list(n_tx = n_tx, ex_starts = ex_starts, ex_ends = ex_ends,
         eq_lens = eq_lens, eq_vals = eq_vals, first_tx = first_tx, last_tx = last_tx)
}

#' @noRd
bambu_next_get_csr <- function(annotations, n_threads = 1L) {
    n_tx <- length(annotations)
    if (n_tx == 0L) return(NULL)
    d <- .bambu_next_extract_digest_args(annotations)
    hit <- .Call(
        "bambu_next_csr_cache_lookup",
        d$n_tx, d$ex_starts, d$ex_ends, d$eq_lens, d$eq_vals, d$first_tx, d$last_tx,
        PACKAGE = "bambu"
    )
    if (!is.null(hit)) return(hit)
    annotations <- bambu_next_attach_csr(annotations, n_threads = n_threads, compute_eq = FALSE)
    d <- .bambu_next_extract_digest_args(annotations)
    .Call(
        "bambu_next_csr_cache_lookup",
        d$n_tx, d$ex_starts, d$ex_ends, d$eq_lens, d$eq_vals, d$first_tx, d$last_tx,
        PACKAGE = "bambu"
    )
}

#' @noRd
bambu_next_attach_csr <- function(annotations, n_threads = 1L, compute_eq = FALSE) {
    n_tx <- length(annotations)
    if (n_tx == 0L) return(annotations)

    d <- .bambu_next_extract_digest_args(annotations)
    if (!compute_eq) {
        hit <- .Call(
            "bambu_next_csr_cache_lookup",
            d$n_tx, d$ex_starts, d$ex_ends, d$eq_lens, d$eq_vals, d$first_tx, d$last_tx,
            PACKAGE = "bambu"
        )
        if (!is.null(hit)) return(annotations)
    }

    ex_counts <- IRanges::elementNROWS(annotations)
    unlisted <- unlist(annotations, use.names = FALSE)
    ex_starts <- d$ex_starts
    ex_ends <- d$ex_ends
    ex_mc <- S4Vectors::mcols(unlisted, use.names = FALSE)
    ex_rank <- if (!is.null(ex_mc$exon_rank)) as.integer(ex_mc$exon_rank) else seq_along(ex_starts)
    ex_end_rank <- if (!is.null(ex_mc$exon_endRank)) as.integer(ex_mc$exon_endRank) else ex_rank

    first_idx <- cumsum(c(1L, ex_counts[-length(ex_counts)]))
    seq_names <- as.character(Seqinfo::seqlevels(annotations))
    seq_lens <- as.integer(Seqinfo::seqlengths(annotations)[seq_names])
    seq_lens[is.na(seq_lens)] <- -1L
    tx_seq_1b <- as.integer(GenomicRanges::seqnames(unlisted))[first_idx]
    st_vec <- as.character(GenomicRanges::strand(unlisted))[first_idx]
    tx_strand_code <- ifelse(st_vec == "+", 1L, ifelse(st_vec == "-", -1L, 0L))

    mc <- S4Vectors::mcols(annotations, use.names = FALSE)
    tx_names <- if (!is.null(mc$TXNAME)) as.character(mc$TXNAME) else as.character(names(annotations))
    tx_gene_names <- if (!is.null(mc$GENEID)) as.character(mc$GENEID) else rep("", n_tx)

    res <- .Call(
        "bambu_next_build_csr_from_flat",
        seq_names, seq_lens, tx_names, tx_gene_names,
        tx_seq_1b, as.integer(tx_strand_code), as.integer(ex_counts),
        as.integer(ex_starts), as.integer(ex_ends),
        as.integer(ex_rank), as.integer(ex_end_rank),
        d$eq_vals, d$eq_lens, isTRUE(compute_eq), as.integer(n_threads),
        PACKAGE = "bambu"
    )
    if (!is.null(res$eq_vals)) {
        part <- IRanges::PartitioningByEnd(cumsum(res$eq_lens), names = NULL)
        S4Vectors::mcols(annotations)$eqClassById <- IRanges::relist(res$eq_vals, part)
        d$eq_lens <- res$eq_lens
        d$eq_vals <- res$eq_vals
    }
    .Call(
        "bambu_next_csr_cache_put",
        res$csr_ptr, d$n_tx, ex_starts, ex_ends, d$eq_lens, d$eq_vals, d$first_tx, d$last_tx,
        PACKAGE = "bambu"
    )
    annotations
}

#' @noRd
bambu_next_prepare_annotations_fast <- function(x, n_threads = 1L) {
    if (is(x, "TxDb")) {
        exonsByTx <- GenomicFeatures::exonsBy(x, by = "tx", use.names = FALSE)
        txNames <- S4Vectors::values(GenomicFeatures::transcripts(x, columns = "tx_name"))$tx_name
        if (any(duplicated(txNames))) {
            message("transcript names are not unique,\n            only one transcript per ID will be kept")
            exonsByTx <- exonsByTx[!duplicated(txNames)]
            txNames <- txNames[!duplicated(txNames)]
        }
        names(exonsByTx) <- txNames
        unlistedExons <- unlist(exonsByTx, use.names = FALSE)
        partitioning <- IRanges::PartitioningByEnd(cumsum(IRanges::elementNROWS(exonsByTx)), names = NULL)
        txIdForReorder <- IRanges::togroup(IRanges::PartitioningByWidth(exonsByTx))
        unlistedExons <- unlistedExons[order(txIdForReorder, unlistedExons$exon_rank)]
        unlistedExons$exon_endRank <- unlist(lapply(IRanges::elementNROWS(exonsByTx), seq, to = 1), use.names = FALSE)
        unlistedExons <- unlistedExons[order(txIdForReorder, IRanges::start(unlistedExons))]
        S4Vectors::mcols(unlistedExons) <- S4Vectors::mcols(unlistedExons)[, c("exon_rank", "exon_endRank")]
        exonsByTx <- IRanges::relist(unlistedExons, partitioning)

        S4Vectors::mcols(exonsByTx) <- suppressMessages(AnnotationDbi::select(
            x, names(exonsByTx), columns = c("TXNAME", "GENEID"), keytype = "TXNAME"
        ))
        S4Vectors::mcols(exonsByTx)$txid <- seq_along(exonsByTx)
        # Compute fast minEqClass via C++ FlatCSR (instead of all-vs-all R findSpliceOverlapsQuick)
        exonsByTx <- bambu_next_attach_csr(exonsByTx, n_threads = n_threads, compute_eq = TRUE)
        return(exonsByTx)
    }

    # GTF file path fast C++ streaming parser
    if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
        stop("Input annotation file not readable. Requires .gtf/.gff format or TxDb object")
    }
    raw_res <- .Call("bambu_next_prepare_gtf_cpp", path.expand(x), as.integer(n_threads), PACKAGE = "bambu")

    # Assemble CompressedGRangesList with exact Bambu S4 slots and attribute order
    strand_map <- ifelse(raw_res$ex_strand_code == 1L, "+",
                  ifelse(raw_res$ex_strand_code == -1L, "-", "*"))
    seq_factor <- factor(raw_res$seq_names[raw_res$ex_seq_idx], levels = raw_res$seq_names)
    st_factor <- factor(strand_map, levels = c("+", "-", "*"))
    gr <- GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(seq_factor),
        ranges = IRanges::IRanges(start = raw_res$ex_start, end = raw_res$ex_end),
        strand = S4Vectors::Rle(st_factor),
        exon_rank = raw_res$ex_rank,
        exon_endRank = raw_res$ex_end_rank
    )
    part_tx <- IRanges::PartitioningByEnd(cumsum(raw_res$tx_ex_counts), names = raw_res$tx_names)
    grlist <- IRanges::relist(gr, part_tx)

    part_eq <- IRanges::PartitioningByEnd(cumsum(raw_res$eq_lens), names = NULL)
    eq_list <- IRanges::relist(raw_res$eq_vals, part_eq)

    mc <- S4Vectors::DataFrame(
        TXNAME = raw_res$tx_names,
        GENEID = raw_res$tx_gene_names,
        txid = seq_along(raw_res$tx_names),
        eqClassById = eq_list
    )
    S4Vectors::mcols(grlist) <- mc
    nt <- length(raw_res$tx_names)
    .Call(
        "bambu_next_csr_cache_put",
        raw_res$csr_ptr, nt, raw_res$ex_start, raw_res$ex_end,
        raw_res$eq_lens, raw_res$eq_vals,
        if (nt > 0L) raw_res$tx_names[1L] else "",
        if (nt > 0L) raw_res$tx_names[nt] else "",
        PACKAGE = "bambu"
    )
    grlist
}
