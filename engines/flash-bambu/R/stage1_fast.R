# Stage 1 C++ Pipeline Glue (Nodes 2 & 3 -> ReadClass SummarizedExperiment before scoring)

#' @noRd
.bambu_next_score_xgb_bundle <- function(bundle, booster) {
    n_j <- length(bundle$junc_to_site)
    if (n_j == 0L) return(numeric(0))
    n_sites <- max(bundle$junc_to_site)
    site_preds <- rep(NA_real_, n_sites)
    if (nrow(bundle$mat) > 0L && !is.null(booster)) {
        p <- xgboost:::predict.xgb.Booster(booster, bundle$mat)
        site_preds[bundle$my_set_all] <- as.numeric(p)
    }
    site_preds[bundle$junc_to_site]
}

#' @noRd
.bambu_next_build_se_from_bsd <- function(
    bsd_ptr,
    csr_ptr,
    stranded,
    verbose,
    run_name,
    ann_seqnames,
    warnings_vec = character(0),
    trackReads = FALSE,
    fallback_ref_seq_names = character(0)
) {
    prep_ec <- .Call("bambu_next_prepare_junc_ec_cpp", bsd_ptr, csr_ptr, isTRUE(stranded), PACKAGE = "bambu")

    can_run_hc <- FALSE
    p_ss <- p_se <- p_es <- p_ee <- numeric(0)
    m_ss <- m_se <- m_es <- m_ee <- numeric(0)

    if (!isTRUE(prep_ec$empty)) {
        if (isTRUE(prep_ec$need_fit)) {
            tb <- prep_ec$train_bundles
            junctionModel <- list(
                spliceSitePredictionStart.start = if (length(tb$Start.start$train_idx) > 0L) {
                    fitXGBoostModel(labels.train = tb$Start.start$train_labels,
                                    data.train = tb$Start.start$mat[tb$Start.start$train_idx, , drop = FALSE],
                                    show.cv = verbose, maxSize.cv = 10000)
                } else NULL,
                spliceSitePredictionStart.end = if (length(tb$Start.end$train_idx) > 0L) {
                    fitXGBoostModel(labels.train = tb$Start.end$train_labels,
                                    data.train = tb$Start.end$mat[tb$Start.end$train_idx, , drop = FALSE],
                                    show.cv = verbose, maxSize.cv = 10000)
                } else NULL,
                spliceSitePredictionEnd.start = if (length(tb$End.start$train_idx) > 0L) {
                    fitXGBoostModel(labels.train = tb$End.start$train_labels,
                                    data.train = tb$End.start$mat[tb$End.start$train_idx, , drop = FALSE],
                                    show.cv = verbose, maxSize.cv = 10000)
                } else NULL,
                spliceSitePredictionEnd.end = if (length(tb$End.end$train_idx) > 0L) {
                    fitXGBoostModel(labels.train = tb$End.end$train_labels,
                                    data.train = tb$End.end$mat[tb$End.end$train_idx, , drop = FALSE],
                                    show.cv = verbose, maxSize.cv = 10000)
                } else NULL
            )
        } else {
            junctionModel <- standardJunctionModels_temp
            if (verbose) message("Junction correction with not enough data, precalculated model is used")
        }

        can_run_hc <- isTRUE(prep_ec$can_run_hc)
        if (can_run_hc) {
            hp <- prep_ec$hc_plus_bundles
            hm <- prep_ec$hc_minus_bundles
            p_ss <- .bambu_next_score_xgb_bundle(hp$Start.start, junctionModel$spliceSitePredictionStart.start)
            p_se <- .bambu_next_score_xgb_bundle(hp$Start.end,   junctionModel$spliceSitePredictionStart.end)
            p_es <- .bambu_next_score_xgb_bundle(hp$End.start,   junctionModel$spliceSitePredictionEnd.start)
            p_ee <- .bambu_next_score_xgb_bundle(hp$End.end,     junctionModel$spliceSitePredictionEnd.end)

            m_ss <- .bambu_next_score_xgb_bundle(hm$Start.start, junctionModel$spliceSitePredictionStart.start)
            m_se <- .bambu_next_score_xgb_bundle(hm$Start.end,   junctionModel$spliceSitePredictionStart.end)
            m_es <- .bambu_next_score_xgb_bundle(hm$End.start,   junctionModel$spliceSitePredictionEnd.start)
            m_ee <- .bambu_next_score_xgb_bundle(hm$End.end,     junctionModel$spliceSitePredictionEnd.end)
        } else {
            warning("no junction correction as no high confidence reference junctions found")
        }
    }

    rc_raw <- .Call(
        "bambu_next_build_read_classes_cpp",
        bsd_ptr, csr_ptr, isTRUE(stranded), can_run_hc,
        p_ss, p_se, p_es, p_ee,
        m_ss, m_se, m_es, m_ee,
        PACKAGE = "bambu"
    )

    meta <- rc_raw$rc_meta
    exons <- rc_raw$exons
    reads <- rc_raw$reads
    seq_names <- rc_raw$seq_names
    ref_seq_names <- if (length(rc_raw$ref_seq_names) > 0L) rc_raw$ref_seq_names else fallback_ref_seq_names

    st_levels <- c("+", "-", "*")
    if (length(meta$rc_names) == 0L) {
        exonsByRC <- GenomicRanges::GRangesList()
        GenomeInfoDb::seqlevels(exonsByRC) <- ref_seq_names
        counts <- matrix(numeric(0), nrow = 0L, ncol = 1L, dimnames = list(NULL, run_name))
        se <- SummarizedExperiment::SummarizedExperiment(
            assays = S4Vectors::SimpleList(counts = counts),
            rowRanges = exonsByRC,
            colData = S4Vectors::DataFrame(name = run_name, row.names = run_name)
        )
        S4Vectors::metadata(se)$warnings <- warnings_vec
        return(se)
    }

    ex_st_char <- ifelse(exons$ex_strand_code == 1L, "+", ifelse(exons$ex_strand_code == -1L, "-", "*"))
    gr <- GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(factor(seq_names[exons$ex_seq_idx], levels = seq_names)),
        ranges = IRanges::IRanges(start = exons$ex_start, end = exons$ex_end),
        strand = S4Vectors::Rle(factor(ex_st_char, levels = st_levels)),
        exon_rank = exons$ex_rank,
        exon_endRank = exons$ex_end_rank
    )
    part_ex <- IRanges::PartitioningByEnd(cumsum(exons$rc_ex_counts), names = meta$rc_names)
    exonsByRC <- IRanges::relist(gr, part_ex)

    part_rid <- IRanges::PartitioningByEnd(cumsum(reads$rc_rid_counts), names = NULL)
    read_ids_list <- as.list(IRanges::relist(reads$flat_read_ids, part_rid))

    conf_labels <- c("highConfidenceJunctionReads", "lowConfidenceJunctionReads",
                     "unsplicedWithin", "unsplicedNew")
    conf_vec <- conf_labels[meta$conf_code]
    has_spliced <- any(meta$conf_code <= 2L)
    has_unspliced <- any(meta$conf_code >= 3L)
    if (!has_spliced) {
        conf_col <- conf_vec
    } else if (!has_unspliced) {
        conf_col <- factor(conf_vec, levels = c("highConfidenceJunctionReads", "lowConfidenceJunctionReads"))
    } else {
        conf_col <- factor(conf_vec, levels = unique(c("highConfidenceJunctionReads", "lowConfidenceJunctionReads", conf_vec)))
    }

    rc_st_char <- ifelse(meta$strand_code == 1L, "+", ifelse(meta$strand_code == -1L, "-", "*"))
    chr_rc_levels <- seq_names
    if (isTRUE(has_unspliced) && isTRUE(rc_raw$all_bam_sq_in_ref) && length(ann_seqnames) > 0L) {
        chr_rc_levels <- unique(c(seq_names, ann_seqnames))
    }
    rc_df <- S4Vectors::DataFrame(
        chr.rc = factor(seq_names[meta$chr_rc_idx], levels = chr_rc_levels),
        strand.rc = factor(rc_st_char, levels = st_levels),
        startSD = meta$startSD,
        endSD = meta$endSD,
        readCount.posStrand = meta$readCount_posStrand,
        intronStarts = meta$intronStarts,
        intronEnds = meta$intronEnds,
        confidenceType = conf_col,
        readCount = meta$readCount
    )
    rc_df$readIds <- read_ids_list
    S4Vectors::mcols(exonsByRC) <- rc_df
    GenomeInfoDb::seqlevels(exonsByRC) <- ref_seq_names

    counts <- matrix(meta$readCount, dimnames = list(meta$rc_names, run_name))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = counts),
        rowRanges = exonsByRC,
        colData = S4Vectors::DataFrame(name = run_name, row.names = run_name)
    )
    S4Vectors::metadata(se)$warnings <- warnings_vec
    if (isTRUE(trackReads) && !is.null(rc_raw$tracked)) {
        S4Vectors::metadata(se)$readNames <- rc_raw$tracked$readNames
        S4Vectors::metadata(se)$readId <- rc_raw$tracked$readId
    }
    se
}

#' @noRd
bambu_next_construct_rc_se_cpp <- function(
    bam_file,
    genomeSequence,
    annotations,
    stranded = FALSE,
    verbose = FALSE,
    trackReads = FALSE,
    n_threads = 4L,
    run_name = NULL,
    lowMemory = FALSE
) {
    bambu_next_apply_thread_guardrails()
    if (is(bam_file, "BamFileList")) {
        if (is.null(run_name) && !is.null(names(bam_file)) && nzchar(names(bam_file)[1L])) {
            run_name <- names(bam_file)[1L]
        }
        bam_file <- bam_file[[1L]]
    }
    if (is(bam_file, "BamFile")) {
        bam_path <- as.character(bam_file$path)[1L]
        bai_raw <- tryCatch(as.character(bam_file$index)[1L], error = function(e) "")
        bai_path <- if (length(bai_raw) == 1L && !is.na(bai_raw) && nzchar(bai_raw) && file.exists(bai_raw)) bai_raw else NULL
        if (is.null(run_name) || is.na(run_name) || !nzchar(run_name)) {
            run_name <- tools::file_path_sans_ext(basename(bam_path))
        }
    } else {
        bam_path <- path.expand(as.character(bam_file)[1L])
        bai_path <- NULL
        if (is.null(run_name) || is.na(run_name) || !nzchar(run_name)) {
            run_name <- tools::file_path_sans_ext(basename(bam_path))
        }
    }

    is_fa <- is(genomeSequence, "FaFile")
    fa_path <- if (is_fa) as.character(genomeSequence$path)[1L] else NULL
    fai_raw <- if (is_fa) tryCatch(as.character(genomeSequence$index)[1L], error = function(e) "") else ""
    fai_path <- if (is_fa && length(fai_raw) == 1L && !is.na(fai_raw) && nzchar(fai_raw) && file.exists(fai_raw)) fai_raw else NULL

    csr_ptr <- if (length(annotations) > 0L) bambu_next_get_csr(annotations, n_threads = n_threads) else NULL

    g_seqnames <- as.character(Seqinfo::seqlevels(genomeSequence))
    g_seqlengths <- as.integer(Seqinfo::seqlengths(genomeSequence)[g_seqnames])
    ann_seqnames <- as.character(Seqinfo::seqlevels(annotations))

    res1 <- .Call(
        "bambu_next_extract_bam_and_motifs_cpp",
        bam_path, bai_path, run_name,
        g_seqnames, g_seqlengths, ann_seqnames,
        fa_path, fai_path, isTRUE(trackReads), as.integer(n_threads),
        PACKAGE = "bambu"
    )
    bsd_ptr <- res1$bsd_ptr

    if (res1$n_spliced == 0L && res1$n_unspliced == 0L) {
        stop("No reads left after filtering.")
    }

    # In-memory BSgenome / DNAStringSet fallback for motifs
    if (!is_fa && res1$n_unique_junc > 0L) {
        uj_info <- .Call("bambu_next_export_unique_junctions_cpp", bsd_ptr, PACKAGE = "bambu")
        used_seq <- uj_info$seq_names[sort(unique(uj_info$seq_idx))]
        seq_chars <- vapply(used_seq, function(sn) {
            as.character(BSgenome::getSeq(genomeSequence, sn))
        }, character(1L))
        .Call("bambu_next_fill_motifs_from_r_cpp", bsd_ptr, used_seq, unname(seq_chars), PACKAGE = "bambu")
    }

    warnings_vec <- as.character(res1$warnings)
    if (verbose && length(warnings_vec) > 0L) {
        warning(paste(warnings_vec, collapse = "\n"))
    }
    if (res1$n_unique_junc == 0L) {
        wtxt <- paste0(
            "No aligned spliced reads detected!",
            "Bambu expects spliced reads. If this is intended, ",
            "see Documentation on how to handle single-exon ",
            "transcripts"
        )
        warnings_vec <- c(warnings_vec, wtxt)
        if (verbose) warning(wtxt)
    }

    if (!isTRUE(lowMemory)) {
        return(.bambu_next_build_se_from_bsd(
            bsd_ptr = bsd_ptr,
            csr_ptr = csr_ptr,
            stranded = stranded,
            verbose = verbose,
            run_name = run_name,
            ann_seqnames = ann_seqnames,
            warnings_vec = warnings_vec,
            trackReads = trackReads
        ))
    }

    # lowMemory = TRUE: construct junction tables & ReadClasses per chromosome cleanly matching native lowMemoryConstructReadClasses
    active_seqs <- .Call("bambu_next_get_bsd_active_seqs_cpp", bsd_ptr, PACKAGE = "bambu")
    se_list <- lapply(active_seqs, function(sid) {
        sub_ptr <- .Call("bambu_next_subset_bsd_seq_cpp", bsd_ptr, sid, PACKAGE = "bambu")
        se_chr <- .bambu_next_build_se_from_bsd(
            bsd_ptr = sub_ptr,
            csr_ptr = csr_ptr,
            stranded = stranded,
            verbose = verbose,
            run_name = run_name,
            ann_seqnames = ann_seqnames,
            warnings_vec = character(0),
            trackReads = FALSE
        )
        if (nrow(se_chr) == 0L) return(NULL)
        se_chr
    })
    se_list <- se_list[!vapply(se_list, is.null, logical(1L))]
    if (length(se_list) == 0L) {
        return(.bambu_next_build_se_from_bsd(
            bsd_ptr = bsd_ptr,
            csr_ptr = csr_ptr,
            stranded = stranded,
            verbose = verbose,
            run_name = run_name,
            ann_seqnames = ann_seqnames,
            warnings_vec = warnings_vec,
            trackReads = trackReads
        ))
    }
    se <- do.call("rbind", se_list)
    if (nrow(se) > 0L) {
        rownames(se) <- paste("rc", seq_len(nrow(se)), sep = ".")
        if (is.character(SummarizedExperiment::rowData(se)$confidenceType) &&
            any(vapply(se_list, function(x) is.factor(SummarizedExperiment::rowData(x)$confidenceType), logical(1L)))) {
            ct_lvl <- unique(unlist(lapply(se_list, function(x) {
                v <- SummarizedExperiment::rowData(x)$confidenceType
                if (is.factor(v)) levels(v) else unique(v)
            })))
            SummarizedExperiment::rowData(se)$confidenceType <- factor(SummarizedExperiment::rowData(se)$confidenceType, levels = ct_lvl)
        }
    }
    S4Vectors::metadata(se)$warnings <- warnings_vec
    if (isTRUE(trackReads)) {
        tr <- .Call("bambu_next_get_bsd_tracked_cpp", bsd_ptr, PACKAGE = "bambu")
        if (!is.null(tr)) {
            S4Vectors::metadata(se)$readNames <- tr$readNames
            S4Vectors::metadata(se)$readId <- tr$readId
        }
    }
    se
}
