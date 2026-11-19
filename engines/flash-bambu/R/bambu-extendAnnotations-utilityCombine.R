#' Combine transcript candidates across samples
#' @param readClassList readClassList
#' @param stranded stranded
#' @param min.readCount minimum number of read counts per transcript
#' @param min.readFractionByGene minimum fraction of transcript usage
#' @param bpParameters biocParallel parameters
#' @param verbose verbose
#' @importFrom GenomicRanges GRanges
#' @importFrom SummarizedExperiment rbind
#' @importFrom dplyr as_tibble mutate_if select mutate %>% filter
#' @importFrom tidyr separate
#' @noRd
.bambu_cast_int_rc <- function(rc_num) {
    if (length(rc_num) > 0L && all(!is.na(rc_num)) &&
        all(floor(rc_num) == rc_num) && max(rc_num) <= 2147483647) {
        return(as.integer(rc_num))
    }
    rc_num
}

.bambu_extract_combine_inputs <- function(readClassList, need_spliced = TRUE, need_unspliced = TRUE) {
    n_sample <- length(readClassList)
    spliced_inputs <- if (need_spliced) vector("list", n_sample) else NULL
    u_seq_lvls <- character(0)
    u_lists <- if (need_unspliced) list(
        seq_names = vector("list", n_sample), st_code = vector("list", n_sample),
        start = vector("list", n_sample), end = vector("list", n_sample),
        rc = vector("list", n_sample), grp = vector("list", n_sample),
        txs = vector("list", n_sample), txsnf = vector("list", n_sample)
    ) else NULL

    for (s in seq_len(n_sample)) {
        rc_obj <- readClassList[[s]]
        if (is.character(rc_obj)) rc_obj <- readRDS(file = rc_obj)
        rr <- SummarizedExperiment::rowRanges(rc_obj)
        rd <- SummarizedExperiment::rowData(rc_obj)
        u_seq_lvls <- unique(c(u_seq_lvls, as.character(Seqinfo::seqlevels(rr))))
        n_rc <- nrow(rc_obj)

        if (n_rc == 0L) {
            if (need_spliced) {
                spliced_inputs[[s]] <- list(
                    ex_counts = integer(0), ex_start = integer(0), ex_end = integer(0),
                    chr = character(0), strand = character(0),
                    intronStarts = character(0), intronEnds = character(0),
                    is_hc = logical(0), readCount = numeric(0), geneReadProp = numeric(0),
                    txScore = numeric(0), txScore_noFit = numeric(0), numExons = integer(0)
                )
            }
            next
        }

        ex_cnts <- as.integer(IRanges::elementNROWS(rr))
        unlisted <- unlist(rr, use.names = FALSE)
        ex_starts <- as.integer(IRanges::start(unlisted))
        ex_ends <- as.integer(IRanges::end(unlisted))
        ct <- as.character(rd$confidenceType)
        rc_num <- as.numeric(rd$readCount)
        grp_num <- as.numeric(rd$geneReadProp)
        txs_num <- as.numeric(rd$txScore)
        txsnf_num <- as.numeric(rd$txScore.noFit)

        if (need_spliced) {
            spliced_inputs[[s]] <- list(
                ex_counts = ex_cnts, ex_start = ex_starts, ex_end = ex_ends,
                chr = as.character(rd$chr.rc), strand = as.character(rd$strand.rc),
                intronStarts = as.character(rd$intronStarts),
                intronEnds = as.character(rd$intronEnds),
                is_hc = (ct == "highConfidenceJunctionReads"),
                readCount = rc_num, geneReadProp = grp_num,
                txScore = txs_num, txScore_noFit = txsnf_num,
                numExons = as.integer(rd$numExons)
            )
        }
        if (need_unspliced) {
            eq_flag <- as.logical(rd$equal)
            eq_flag[is.na(eq_flag)] <- FALSE
            idx_u <- which(ct == "unsplicedNew" & rc_num > 1 & !eq_flag)
            if (length(idx_u) > 0L) {
                first_pos <- cumsum(c(1L, ex_cnts[-n_rc]))[idx_u]
                st_str <- as.character(GenomicRanges::strand(unlisted))[first_pos]
                u_lists$seq_names[[s]] <- as.character(GenomicRanges::seqnames(unlisted))[first_pos]
                u_lists$st_code[[s]] <- ifelse(st_str == "+", 1L, ifelse(st_str == "-", -1L, 0L))
                u_lists$start[[s]] <- ex_starts[first_pos]
                u_lists$end[[s]] <- ex_ends[first_pos]
                u_lists$rc[[s]] <- rc_num[idx_u]
                u_lists$grp[[s]] <- grp_num[idx_u]
                u_lists$txs[[s]] <- txs_num[idx_u]
                u_lists$txsnf[[s]] <- txsnf_num[idx_u]
            }
        }
    }
    list(spliced = spliced_inputs, unspliced = u_lists, seq_lvls = u_seq_lvls)
}

.bambu_run_kway_spliced <- function(spliced_inputs, seq_lvls, min.readCount,
        min.readFractionByGene, min.txScore.multiExon, min.txScore.singleExon) {
    sp_res <- .Call(
        "bambu_next_kway_merge_spliced_cpp",
        spliced_inputs,
        as.numeric(min.readCount), as.numeric(min.readFractionByGene),
        as.numeric(min.txScore.multiExon), as.numeric(min.txScore.singleExon),
        PACKAGE = "bambu"
    )
    chr_lvls <- if (length(seq_lvls) > 0L) seq_lvls else unique(sp_res$chr)
    data.table(
        start = sp_res$start,
        end = sp_res$end,
        readCount = .bambu_cast_int_rc(sp_res$readCount),
        intronStarts = sp_res$intronStarts,
        intronEnds = sp_res$intronEnds,
        chr = factor(sp_res$chr, levels = chr_lvls),
        strand = factor(sp_res$strand, levels = c("+", "-", "*")),
        maxTxScore = sp_res$maxTxScore,
        maxTxScore.noFit = sp_res$maxTxScore_noFit,
        NSampleReadCount = sp_res$NSampleReadCount,
        NSampleReadProp = sp_res$NSampleReadProp,
        NSampleTxScore = sp_res$NSampleTxScore
    )
}

.bambu_run_merge_unspliced <- function(u_lists, seq_lvls, stranded, min.readCount,
        min.readFractionByGene, min.txScore.singleExon) {
    u_seq_all <- unlist(u_lists$seq_names, use.names = FALSE)
    if (length(u_seq_all) == 0L) {
        return(data.table(
            chr = factor(character(0), levels = seq_lvls),
            strand = factor(character(0), levels = c("+", "-", "*")),
            start = integer(0), end = integer(0), readCount = integer(0),
            maxTxScore = numeric(0), maxTxScore.noFit = numeric(0),
            NSampleReadCount = integer(0), NSampleReadProp = integer(0), NSampleTxScore = integer(0)
        ))
    }
    un_res <- .Call(
        "bambu_next_merge_unspliced_cpp",
        seq_lvls,
        as.integer(match(u_seq_all, seq_lvls)),
        as.integer(unlist(u_lists$st_code, use.names = FALSE)),
        as.integer(unlist(u_lists$start, use.names = FALSE)),
        as.integer(unlist(u_lists$end, use.names = FALSE)),
        as.numeric(unlist(u_lists$rc, use.names = FALSE)),
        as.numeric(unlist(u_lists$grp, use.names = FALSE)),
        as.numeric(unlist(u_lists$txs, use.names = FALSE)),
        as.numeric(unlist(u_lists$txsnf, use.names = FALSE)),
        isTRUE(stranded),
        as.numeric(min.readCount), as.numeric(min.readFractionByGene),
        as.numeric(min.txScore.singleExon),
        PACKAGE = "bambu"
    )
    st_map <- ifelse(un_res$strand_code == 1L, "+", ifelse(un_res$strand_code == -1L, "-", "*"))
    data.table(
        chr = factor(seq_lvls[un_res$chr_idx], levels = seq_lvls),
        strand = factor(st_map, levels = c("+", "-", "*")),
        start = un_res$start,
        end = un_res$end,
        readCount = .bambu_cast_int_rc(un_res$readCount),
        maxTxScore = un_res$maxTxScore,
        maxTxScore.noFit = un_res$maxTxScore_noFit,
        NSampleReadCount = un_res$NSampleReadCount,
        NSampleReadProp = un_res$NSampleReadProp,
        NSampleTxScore = un_res$NSampleTxScore
    )
}

isore.combineTranscriptCandidates <- function(readClassList,
    stranded, ## stranded used for unspliced reduce  
    min.readCount , min.readFractionByGene,
    min.txScore.multiExon, min.txScore.singleExon, bpParameters ,verbose){
    start.ptm <- proc.time()
    ext <- .bambu_extract_combine_inputs(readClassList, TRUE, TRUE)
    combinedSplicedTranscripts <- .bambu_run_kway_spliced(
        ext$spliced, ext$seq_lvls, min.readCount, min.readFractionByGene,
        min.txScore.multiExon, min.txScore.singleExon)
    combinedSplicedTranscripts[, confidenceType := "highConfidenceJunctionReads"]
    combinedUnsplicedTranscripts <- .bambu_run_merge_unspliced(
        ext$unspliced, ext$seq_lvls, stranded, min.readCount,
        min.readFractionByGene, min.txScore.singleExon)
    combinedUnsplicedTranscripts[, confidenceType := "unsplicedNew"]
    combinedTranscripts <- as_tibble(rbindlist(list(combinedSplicedTranscripts,
        combinedUnsplicedTranscripts), fill = TRUE))
    end.ptm <- proc.time()
    if (verbose) message("combined spliced/unspliced candidates in ",
        round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    return(combinedTranscripts)
}


#' combine spliced transcript models
#' @noRd
combineSplicedTranscriptModels <- function(readClassList, bpParameters, 
        min.readCount, min.readFractionByGene, min.txScore.multiExon, 
        min.txScore.singleExon, verbose){
    options(scipen = 999)
    start.ptm <- proc.time()
    ext <- .bambu_extract_combine_inputs(readClassList, TRUE, FALSE)
    combinedFeatureTibble <- .bambu_run_kway_spliced(
        ext$spliced, ext$seq_lvls, min.readCount, min.readFractionByGene,
        min.txScore.multiExon, min.txScore.singleExon)
    end.ptm <- proc.time()
    if (verbose) message("combing spliced feature tibble objects across all ",
        "samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(combinedFeatureTibble)
}

#' Sequentially combine feature tibbles 
#' @noRd
sequentialCombineFeatureTibble <- function(readClassList,
        indexList,intraGroup,min.readCount,min.readFractionByGene,
        min.txScore.multiExon, min.txScore.singleExon){
    combinedFeatureTibble <- NULL
    for (s in seq_along(readClassList)){
        combinedListNew <- readClassList[[s]]
        if(intraGroup){
            combinedListNew <- 
                extractFeaturesFromReadClassSE(readClassSe = combinedListNew,
                    sample_id = indexList[s], min.readCount = min.readCount,
                    min.readFractionByGene = min.readFractionByGene,
                    min.txScore.multiExon = min.txScore.multiExon,
                    min.txScore.singleExon = min.txScore.singleExon)
        }
        combinedFeatureTibble <- combineFeatureTibble(combinedFeatureTibble,
            combinedListNew, index = indexList[s], intraGroup)
    }
    return(combinedFeatureTibble)
}



#' @noRd 
updateStartEndReadCount <- function(combinedFeatureTibble){
    setDT(combinedFeatureTibble)
    combinedFeatureTibble[, rowID := .I]
    
    colNames <- colnames(combinedFeatureTibble)
    readCountCols <- sort(colNames[grep("^readCount", colNames)]) # to make sure it's ordered by sample name
    startCols <- sort(colNames[grep("^start", colNames)])
    endCols <- sort(colNames[grep("^end", colNames)])
    
    startEndDt <- combinedFeatureTibble[, 
        .(start = readCountWeightedMedian(.SD,x,y),
        end = readCountWeightedMedian(.SD,z,y),
        readCount = sum(.SD[,y], na.rm = TRUE)),
        by = rowID,  env = I(list(x = startCols, y = readCountCols,z = endCols))]

    combinedFeatureTibble <- startEndDt[combinedFeatureTibble[,.(intronStarts, intronEnds, chr, strand, maxTxScore, 
                                                                 maxTxScore.noFit, NSampleReadCount, NSampleReadProp, 
                                                                 NSampleTxScore, rowID)], on = "rowID"]
    combinedFeatureTibble[, rowID := NULL]
    return(combinedFeatureTibble)
}

#' Function to get median value without interpolation using certain column names
#' @noRd
readCountWeightedMedian <- function(dt, valuevar, timesvar){
    sortVector <- rep(na.omit(unlist(dt[,..valuevar])), 
                times = as.integer(na.omit(unlist(dt[,..timesvar]))))
    return(min(sortVector[sortVector>=quantile(sortVector, probs = 0.5)]))
}


#' Function to combine featureTibble and create the NSample variables 
#' @noRd
combineFeatureTibble <- function(combinedFeatureTibble,
        featureTibbleSummarised, index=1, intraGroup = TRUE){ 
    if (is.null(combinedFeatureTibble)) { 
        combinedTable <- featureTibbleSummarised %>% 
            select(intronStarts, intronEnds, chr, strand, maxTxScore, 
            maxTxScore.noFit, NSampleReadCount, NSampleReadProp,NSampleTxScore, 
            starts_with('start'), starts_with('end'), starts_with('readCount'))
    } else { 
        combinedTable = full_join(combinedFeatureTibble, 
            featureTibbleSummarised, by = c('intronStarts', 'intronEnds', 'chr',
            'strand'), suffix=c('.combined','.new')) %>% 
            mutate(NSampleReadCount=pmax0NA(NSampleReadCount.combined) + 
                        pmax0NA(NSampleReadCount.new), 
                    NSampleReadProp = pmax0NA(NSampleReadProp.combined) + 
                        pmax0NA(NSampleReadProp.new), 
                    NSampleTxScore = pmax0NA(NSampleTxScore.combined) + 
                        pmax0NA(NSampleTxScore.new),
                    maxTxScore = pmax(maxTxScore.combined, 
                        maxTxScore.new, na.rm = TRUE),
                    maxTxScore.noFit = pmax(maxTxScore.noFit.combined, 
                        maxTxScore.noFit.new, na.rm = TRUE)) %>% 
            select(intronStarts, intronEnds, chr, strand,
            NSampleReadCount, NSampleReadProp, NSampleTxScore, maxTxScore, 
            maxTxScore.noFit, starts_with('start'), starts_with('end'), 
            starts_with('readCount')) 
    } 
    if(intraGroup) 
        combinedTable <- 
            rename_with(combinedTable, ~gsub('^(end|start|readCount)$',
                                            paste0('\\1\\.',index), .x)) 
    return(combinedTable) 
}
#' pmax replace NAs with 0
#' @noRd
pmax0NA <- function(vec){
    vec[is.na(vec)] <- 0
    return(pmax(vec))
}
#' pmin replace NAs with 0
#' @noRd
pmin0NA <- function(vec){
    vec[is.na(vec)] <- 0
    return(pmin(vec))
}
#' extract important features from readClassSe object for each sample
#' @noRd
extractFeaturesFromReadClassSE <- function(readClassSe, sample_id,
        min.readCount, min.readFractionByGene, 
        min.txScore.multiExon, min.txScore.singleExon){
    if (is.character(readClassSe)) 
        readClassSe <- readRDS(file = readClassSe)
    dimNames <- list(rownames(readClassSe), colnames(readClassSe))
    rowRangesSe <- rowRanges(readClassSe)
    rowData <- as_tibble(rowData(readClassSe)) %>% 
        mutate(start = unname(min(start(rowRangesSe))), 
                end= unname(max(end(rowRangesSe))))
    group_var <- c("intronStarts", "intronEnds", "chr", "strand")
    sum_var <- c("start","end","NSampleReadCount", "maxTxScore", 
                "maxTxScore.noFit", "readCount","NSampleReadProp",
                "NSampleTxScore")
    featureTibble <- rowData %>% 
        dplyr::select(chr = chr.rc, start, end, strand = strand.rc, 
            intronStarts, intronEnds, confidenceType, readCount, geneReadProp, 
            txScore, txScore.noFit, numExons) %>%
        filter(readCount >= 1, # only use readCount>1 and highconfidence reads
            confidenceType == "highConfidenceJunctionReads") %>% 
        mutate(NSampleReadCount = (readCount >= min.readCount), 
            # number of samples passed read count criteria
            NSampleReadProp = (geneReadProp >= min.readFractionByGene),
            NSampleTxScore = ((txScore > min.txScore.multiExon & numExons >= 2) |
            (txScore > min.txScore.singleExon & numExons == 1)), 
            maxTxScore = txScore, maxTxScore.noFit = txScore.noFit) %>%
        select(all_of(c(group_var, sum_var))) 
    return(featureTibble)
}


#' combine unspliced transcript models
#' @importFrom tidyr separate
#' @importFrom dplyr %>% select mutate 
#' @importFrom summarizedExperiment rowRanges
#' @importFrom biocParallel bplapply
#' @noRd
combineUnsplicedTranscriptModels <- 
    function(readClassList,  bpParameters, stranded, min.readCount, 
            min.readFractionByGene, min.txScore.multiExon,
            min.txScore.singleExon, verbose){
        start.ptm <- proc.time()
        ext <- .bambu_extract_combine_inputs(readClassList, FALSE, TRUE)
        combinedUnsplicedTibble <- as_tibble(.bambu_run_merge_unspliced(
            ext$unspliced, ext$seq_lvls, stranded, min.readCount,
            min.readFractionByGene, min.txScore.singleExon))
        end.ptm <- proc.time()
        if (verbose) message("combine new unspliced tibble object across all ",
        "samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
        return(combinedUnsplicedTibble)
    }


#' extract new unspliced ranges from readClassSe object for each sample
#' @importFrom dplyr as_tibble
#' @importFrom summarizedExperiment rowData rownames
#' @noRd
extractNewUnsplicedRanges <- function(readClassSe, sample_id){
    if (is.character(readClassSe)) 
        readClassSe <- readRDS(file = readClassSe)
    rowData <- as_tibble(rowData(readClassSe))
    pre_names <- rownames(readClassSe)
    rownames(readClassSe) <- paste0("s",sample_id,"-",pre_names)
    newUnsplicedSe <- 
        readClassSe[which(rowData$confidenceType == "unsplicedNew" & 
                                rowData$readCount > 1 &
                                (!rowData$equal))]
    return(newUnsplicedSe)
}


#' reduce unspliced ranges
#' @importFrom dplyr as_tibble %>% mutate group_by summarise ungroup
#' @noRd
reduceUnsplicedRanges <- function(rangesList, stranded){
    unlistedSe <- do.call("c",rangesList)
    combinedNewUnsplicedSe <- 
        reduce(unlistedSe, with.revmap=TRUE,ignore.strand = !stranded)
    ## map it back to find the corresponding sample and id, update name
    rcNames <- as_tibble(as.data.frame(combinedNewUnsplicedSe$revmap)) %>%
        mutate(row_id = unlistedSe$row_id[value]) %>%
        group_by(group) %>%
        summarise(combinedName = paste(row_id, collapse = "+")) %>%
        ungroup()
    ## at each iteration update the new combined name of unspliced se 
    combinedNewUnsplicedSe$row_id <- rcNames$combinedName
    combinedNewUnsplicedSe$revmap <- NULL
    return(combinedNewUnsplicedSe)
}

#' make unspliced tibble 
#' @importFrom tidyr separate_rows pivot_wider
#' @importFrom dplyr as_tibble rename mutate select %>% group_by left_join
#'              ungroup
#' @noRd
makeUnsplicedTibble <- function(combinedNewUnsplicedSe,newUnsplicedSeList,
        colDataNames,min.readCount, min.readFractionByGene,
        min.txScore.multiExon, min.txScore.singleExon, bpParameters){
        bpParameters$progressbar = FALSE
    newUnsplicedTibble <- as_tibble(combinedNewUnsplicedSe) %>%
        rename(chr = seqnames) %>% select(chr, start, end, strand, row_id) %>%
        separate_rows(row_id, sep = "\\+") 
    rowDataCombined <-
        do.call("rbind",bplapply(newUnsplicedSeList, function(newUnsplicedSe) {
            rr <- rowData(newUnsplicedSe[intersect(rownames(newUnsplicedSe), 
                                        newUnsplicedTibble$row_id)])
            rr <- as_tibble(rr) %>% select(confidenceType,readCount, 
                    geneReadProp, txScore, txScore.noFit) %>%
                mutate(row_id = rownames(rr))
            return(rr)
        } , BPPARAM = bpParameters))
    newUnsplicedTibble <- newUnsplicedTibble %>% 
        left_join(rowDataCombined, by =  "row_id") %>%
        mutate(readCount_tmp = readCount) %>%
        group_by(chr,strand, start, end) %>%
        summarise(readCount = sum(readCount),
                  maxTxScore = weighted.mean(txScore, readCount_tmp),
                  maxTxScore.noFit = weighted.mean(txScore.noFit, readCount_tmp),
                  NSampleReadCount = sum(readCount_tmp >= min.readCount), 
                  NSampleReadProp = sum(geneReadProp >= 
                                          min.readFractionByGene),
                  NSampleTxScore = sum(txScore > min.txScore.singleExon))
    
    return(newUnsplicedTibble)
}
