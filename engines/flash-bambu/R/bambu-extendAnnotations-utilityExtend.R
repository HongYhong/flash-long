#' Extend annotations
#' @inheritParams bambu
#' @noRd
isore.extendAnnotations <- function(combinedTranscripts, annotationGrangesList,
                                    remove.subsetTx = TRUE,
                                    min.sampleNumber = 1, NDR = NULL, min.exonDistance = 35, min.exonOverlap = 10,
                                    min.primarySecondaryDist = 5, min.primarySecondaryDistStartEnd = 5, 
                                    min.readFractionByEqClass = 0, fusionMode = FALSE,
                                    prefix = "Bambu", baselineFDR = 0.1, defaultModels = NULL,
                                    verbose = FALSE, ncore = 1L){
  combinedTranscripts <- filterTranscripts(combinedTranscripts, min.sampleNumber)
  if (nrow(combinedTranscripts) > 0) {
    ref_csr_ptr <- if (length(annotationGrangesList) > 0L) {
        bambu_next_get_csr(annotationGrangesList, n_threads = ncore)
    } else NULL
    group_var <- c("intronStarts","intronEnds","chr","strand","start","end",
                   "confidenceType","readCount", "maxTxScore", "maxTxScore.noFit")
    rowDataTibble <- select(combinedTranscripts,all_of(group_var))
    annotationSeqLevels <- seqlevels(annotationGrangesList)
    rowDataSplicedTibble <- filter(rowDataTibble,
                                   confidenceType == "highConfidenceJunctionReads")
    transcriptRanges <- makeExonsIntronsSpliced(
      rowDataSplicedTibble, annotationSeqLevels, build_introns = FALSE)
    confidenceTypeVec <- rowDataTibble$confidenceType
    if(nrow(rowDataSplicedTibble)>0){
        rowDataFilteredSpliced <- addNewSplicedReadClasses(transcriptRanges,
                                                       rowDataSplicedTibble, annotationGrangesList, 
                                                       min.exonDistance, min.primarySecondaryDist,
                                                       min.primarySecondaryDistStartEnd, verbose,
                                                       ref_csr_ptr = ref_csr_ptr, ncore = ncore)
    } else{ rowDataFilteredSpliced = NULL}
    rowDataFilteredUnspliced <- rowDataTibble[which(confidenceTypeVec == "unsplicedNew"),]
    SEnRng <- addNewUnsplicedReadClasses(rowDataFilteredUnspliced, 
                                         rowDataFilteredSpliced, transcriptRanges$exons, 
                                         annotationGrangesList, min.exonOverlap, verbose,
                                         ref_csr_ptr = ref_csr_ptr)
    rowDataCombined <- SEnRng$rowDataCombined
    exonRangesCombined <- SEnRng$exonRangesCombined
    # assign gene IDs based on exon match
    geneIds <- assignGeneIds(grl = exonRangesCombined,
                             annotations = annotationGrangesList,
                             min.exonOverlap = min.exonOverlap,
                             fusionMode = fusionMode)
    rowDataCombined$GENEID <- geneIds[,1]
    rowDataCombined$novelGene <- geneIds[,2]
    if(fusionMode) rowDataCombined$readClassType[geneIds[,3]] <- 'fusionTranscript'
    
    # ## filter out transcripts
    extendedAnnotationRanges <- filterTranscriptsByAnnotation(
      rowDataCombined, annotationGrangesList, exonRangesCombined, prefix,
      remove.subsetTx, min.readFractionByEqClass, baselineFDR, NDR, defaultModels, verbose,
      ref_csr_ptr = ref_csr_ptr, ncore = ncore)
    return(extendedAnnotationRanges)
  } else {
    message("The current filtering criteria filters out all new read 
            classes, please consider less strigent critria!")
    return(annotationGrangesList)
  }
}

#' filter transcripts by read counts
#' @noRd
filterTranscripts <- function(combinedTranscripts, min.sampleNumber){
  filterSet = NULL
  if (nrow(combinedTranscripts) > 0){
    #filter by read counts
    naTxscore <- is.na(combinedTranscripts$NSampleTxScore)
    if(any(naTxscore)) 
       combinedTranscripts$NSampleTxScore[naTxscore] <- 0
    filterSet <- combinedTranscripts$NSampleReadCount >= min.sampleNumber & (
        combinedTranscripts$NSampleTxScore >= min.sampleNumber) & (
        combinedTranscripts$NSampleReadProp >= min.sampleNumber)
  }
  combinedTranscripts = combinedTranscripts[filterSet,]
  return(combinedTranscripts)
}

#' calculate minimum equivalent classes for extended annotations
#' @importFrom dplyr select as_tibble %>% mutate_at mutate group_by 
#'     ungroup .funs .name_repair vars 
#' @noRd
filterTranscriptsByAnnotation <- function(rowDataCombined, annotationGrangesList,
                                          exonRangesCombined, prefix,  remove.subsetTx, 
                                          min.readFractionByEqClass, baselineFDR = 0.1, 
                                          NDR = NULL, defaultModels = NULL, verbose,
                                          ref_csr_ptr = NULL, ncore = 1L) {
  start.ptm <- proc.time() # (1) based on transcript usage
  
  #calculate relative read count before any filtering
  rowDataCombined <- group_by(rowDataCombined, GENEID) %>% mutate(relReadCount = readCount/sum(readCount))
  if (remove.subsetTx) { # (1) based on compatibility with annotations
    notCompatibleIds <- which(!grepl("compatible", rowDataCombined$readClassType) |
        rowDataCombined$readClassType == "equal:compatible") #keep equal for FDR calculation
    exonRangesCombined <- exonRangesCombined[notCompatibleIds]
    rowDataCombined <- rowDataCombined[notCompatibleIds,]
  }
  #(2) remove transcripts below NDR threshold/identical junctions to annotations
  rowDataCombined = calculateNDROnTranscripts(rowDataCombined, 
                        useTxScore = length(annotationGrangesList)==0)

  if(length(annotationGrangesList)>0){ #only recommend an NDR if its possible to calculate an NDR
      NDR = recommendNDR(rowDataCombined, baselineFDR, NDR, defaultModels, verbose)
  } else {
      if(is.null(NDR)) NDR = 0.5
  }
  filterSet = (rowDataCombined$NDR <= NDR)
  exonRangesCombined <- exonRangesCombined[filterSet]
  rowDataCombined <- rowDataCombined[filterSet,]
  #calculate relative subset read count after filtering (increase speed, subsets are not considered here)
  mcols(exonRangesCombined)$txid <- seq_along(exonRangesCombined)
  if (length(exonRangesCombined) > 0L) {
    S4Vectors::mcols(exonRangesCombined)$eqClassById <- NULL
    exonRangesCombined <- bambu_next_attach_csr(exonRangesCombined, n_threads = ncore, compute_eq = TRUE)
    minEq <- as.list(S4Vectors::mcols(exonRangesCombined)$eqClassById)
    rc_vec <- rowDataCombined$readCount
    rowDataCombined$relSubsetCount <- rc_vec / vapply(minEq, function(x) sum(rc_vec[x]), numeric(1))
  } else {
    rowDataCombined$relSubsetCount <- numeric(0)
  }
  #post extend annotation filters applied here (currently only subset filter)
  if(min.readFractionByEqClass>0 & sum(filterSet)>0) { # filter out subset transcripts based on relative expression
    filterSet <- rowDataCombined$relSubsetCount > min.readFractionByEqClass
    exonRangesCombined <- exonRangesCombined[filterSet]
    rowDataCombined <- rowDataCombined[filterSet,]
  }
  if(sum(filterSet)==0 & length(annotationGrangesList)==0) stop(
    "WARNING - No annotations were provided. Please increase NDR threshold to use novel transcripts")
  if(sum(filterSet)==0) message("WARNING - No novel transcripts meet the given thresholds. Try a higher NDR.")
  # (3) combine novel transcripts with annotations
  extendedAnnotationRanges <- combineWithAnnotations(
    rowDataCombined, exonRangesCombined, 
    annotationGrangesList, prefix,
    ref_csr_ptr = ref_csr_ptr, ncore = ncore)
  metadata(extendedAnnotationRanges)$NDRthreshold = NDR
  end.ptm <- proc.time()
  if (verbose) message("transcript filtering in ",
                       round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  return(extendedAnnotationRanges)
}

#' calculates an expected NDR based on the annotations'
#' @noRd
recommendNDR <- function(combinedTranscripts, baselineFDR = 0.1, NDR = NULL, defaultModels = defaultModels, verbose = FALSE){
    if(verbose) message("-- Predicting annotation completeness to determine NDR threshold --")
    equal = combinedTranscripts$readClassType == "equal:compatible"
    equal[is.na(equal)] = FALSE

    #add envirnment so poly() works
	if(is.null(defaultModels$lmNDR[["terms"]])){
		frm <- defaultModels$lmNDR$call$formula
        defaultModels$lmNDR[["terms"]] <- terms(as.formula(frm))
        attr(defaultModels$lmNDR[["terms"]], ".Environment") <- new.env(parent = parent.env(globalenv()))
	}
    #attr(defaultModels$lmNDR[["terms"]], ".Environment") <- new.env(parent = parent.env(globalenv()))
    baseline = predict(defaultModels$lmNDR, newdata=data.frame(NDR=baselineFDR))
    attr(defaultModels$lmNDR[["terms"]], ".Environment") = c()

    NDRscores = calculateNDR(combinedTranscripts$maxTxScore.noFit, equal)
    score = combinedTranscripts$maxTxScore.noFit
    NDR.rec = predict(lm(NDRscores~poly(score,3,raw=TRUE)), newdata=data.frame(score=baseline))
    NDR.rec = round(NDR.rec,3)
    if(verbose) message("Recommended NDR for baseline FDR of ", baselineFDR, " = ", NDR.rec)
    if (NDR.rec < 0) NDR.rec = 0
    if(NDR.rec > 0.5){
        message("A high NDR threshold is being recommended by Bambu indicating high levels of novel transcripts, ",
        "limiting the performance of the trained model")
        message("We recommend training a new model on similiar but well annotated dataset if available ",
        "(https://github.com/GoekeLab/bambu/tree/master#Training-a-model-on-another-speciesdataset-and-applying-it)",
        ", or alternatively running Bambu with opt.discovery=list(fitReadClassModel=FALSE)")
    }
    
    #if users are using an NDR let them know if the recommended NDR is different
    if(is.null(NDR)) 
    {
        NDR = NDR.rec
        message("Using a novel discovery rate (NDR) of: ", NDR)
    }
    else{
        if(abs(NDR.rec-NDR)>=0.1){
            message(paste0("For your combination of sample and reference annotations we recommend an NDR of ", NDR.rec,
            ". You are currently using an NDR threshold of ", NDR, 
            ". A higher NDR is suited for samples where the reference annotations are poor and more novel transcripts are expected,", 
            "whereas a lower NDR is suited for samples with already high quality annotations"))
        }
    }
    return(NDR)
}


#' Calculate NDR based on transcripts 
#' @noRd
calculateNDROnTranscripts <- function(combinedTranscripts, useTxScore = FALSE){
      # calculate and filter by NDR
    equal = combinedTranscripts$readClassType == "equal:compatible"
    equal[is.na(equal)] = FALSE
    if(sum(equal, na.rm = TRUE)<50 | sum(!equal, na.rm = TRUE)<50 | useTxScore){
          combinedTranscripts$NDR = 1 - combinedTranscripts$maxTxScore
          if(!useTxScore) message("WARNING - Less than 50 TRUE or FALSE read classes ",
            "for NDR precision stabilization.")
          message("NDR will be approximated as: (1 - Transcript Model Prediction Score)")
    } else combinedTranscripts$NDR = calculateNDR(combinedTranscripts$maxTxScore, equal)
    return(combinedTranscripts)
}

#' calculates the minimum NDR for each score 
#' @noRd
calculateNDR = function(score, labels){
    scoreOrder = order(score, decreasing = TRUE) 
    labels = labels[scoreOrder]
    NDR = cumsum(!labels)/(seq_len(length(labels))) #calculate NDR
    NDR = rev(cummin(rev(NDR))) #flatten NDR so its never higher than a lower ranked RC
    return(NDR[order(scoreOrder)]) #return to original order
}

#' generate exon/intron ByReadClass objects
#' @importFrom GenomicRanges makeGRangesListFromFeatureFragments
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom dplyr select distinct 
#' @noRd
makeExonsIntronsSpliced <- function(transcriptsTibble, annotationSeqLevels, build_introns = TRUE){
  if (nrow(transcriptsTibble) == 0L) {
      ex <- GRangesList()
      seqlevels(ex) <- unique(c(seqlevels(ex), annotationSeqLevels))
      intr <- GRangesList()
      seqlevels(intr) <- unique(c(seqlevels(intr), annotationSeqLevels))
      return(list("exons" = ex, "introns" = intr, "flat" = NULL))
  }
  tx_chrs <- as.character(transcriptsTibble$chr)
  chr_init <- if (is.factor(transcriptsTibble$chr)) levels(transcriptsTibble$chr) else unique(tx_chrs)
  all_seq_lvls <- unique(c(chr_init, annotationSeqLevels))
  tx_seq_1b <- match(tx_chrs, all_seq_lvls)
  tx_st_str <- as.character(transcriptsTibble$strand)
  tx_st_code <- ifelse(tx_st_str == "+", 1L, ifelse(tx_st_str == "-", -1L, 0L))

  flat <- .Call(
      "bambu_next_build_exons_from_introns_cpp",
      as.integer(tx_seq_1b),
      as.integer(transcriptsTibble$start),
      as.integer(transcriptsTibble$end),
      as.integer(tx_st_code),
      as.character(transcriptsTibble$intronStarts),
      as.character(transcriptsTibble$intronEnds),
      PACKAGE = "bambu"
  )

  st_map <- ifelse(flat$ex_strand_code == 1L, "+",
            ifelse(flat$ex_strand_code == -1L, "-", "*"))
  ex_gr <- GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(factor(all_seq_lvls[flat$ex_seq_1b], levels = all_seq_lvls)),
      ranges = IRanges::IRanges(start = flat$ex_start, end = flat$ex_end),
      strand = S4Vectors::Rle(factor(st_map, levels = c("+", "-", "*"))),
      exon_rank = as.integer(flat$ex_rank),
      exon_endRank = as.integer(flat$ex_end_rank)
  )
  ex_part <- IRanges::PartitioningByEnd(
      cumsum(flat$tx_ex_counts),
      names = as.character(seq_len(nrow(transcriptsTibble)))
  )
  exonsByReadClass <- IRanges::relist(ex_gr, ex_part)

  intronsByReadClass <- GRangesList()
  if (isTRUE(build_introns) && !all(is.na(transcriptsTibble$intronStarts))) {
      intronsByReadClass <- makeGRangesListFromFeatureFragments(
          seqnames = transcriptsTibble$chr,
          fragmentStarts = transcriptsTibble$intronStarts,
          fragmentEnds = transcriptsTibble$intronEnds,
          strand = transcriptsTibble$strand
      )
      names(intronsByReadClass) <- seq_along(intronsByReadClass)
      seqlevels(intronsByReadClass) <- unique(c(seqlevels(intronsByReadClass), annotationSeqLevels))
  }

  flat$seq_names <- all_seq_lvls
  return(list("exons" = exonsByReadClass, "introns" = intronsByReadClass, "flat" = flat))
}


#' create exonsByReadClass
#' @param seFilteredSpliced a SummarizedExperiment object
#' for filtered spliced reads
#' @param annotationGrangesList annotation GRangesList object
#' @importFrom GenomicRanges makeGRangesListFromFeatureFragments unlist narrow
#'     elementNROWS relist
#' @importFrom GenomeInfoDb seqlevels
#' @noRd
createExonByReadClass <- function(transcriptsTibble, annotationSeqLevels) {
  makeExonsIntronsSpliced(transcriptsTibble, annotationSeqLevels, build_introns = FALSE)$exons
}



#' extended annotations for spliced reads
#' @noRd
addNewSplicedReadClasses <- function(combinedTranscriptRanges, 
                                     rowDataFilteredSpliced, annotationGrangesList, min.exonDistance, 
                                     min.primarySecondaryDist, min.primarySecondaryDistStartEnd, verbose,
                                     ref_csr_ptr = NULL, ncore = 1L){
  start.ptm <- proc.time()
  exonsByReadClass <- combinedTranscriptRanges$exons
  flat <- combinedTranscriptRanges$flat
  if (is.null(flat)) {
      flat <- .bambu_extract_grl_flat(exonsByReadClass, combinedTranscriptRanges$introns)
  }
  if (is.null(ref_csr_ptr) && length(annotationGrangesList) > 0L) {
      ref_csr_ptr <- bambu_next_get_csr(annotationGrangesList, n_threads = ncore)
  }
  res <- .Call(
      "bambu_next_add_new_spliced_cpp",
      ref_csr_ptr,
      as.character(flat$seq_names),
      as.integer(flat$tx_ex_counts),
      as.integer(flat$ex_seq_1b),
      as.integer(flat$ex_start),
      as.integer(flat$ex_end),
      as.integer(flat$ex_strand_code),
      as.integer(flat$ex_rank),
      as.integer(flat$ex_end_rank),
      as.integer(flat$tx_in_counts),
      as.integer(flat$in_start),
      as.integer(flat$in_end),
      as.integer(min.exonDistance),
      as.integer(min.primarySecondaryDist),
      as.integer(min.primarySecondaryDistStartEnd),
      as.integer(ncore),
      PACKAGE = "bambu"
  )
  rowDataFilteredSpliced$TXNAME <- res$TXNAME
  rowDataFilteredSpliced$GENEID <- res$GENEID
  rowDataFilteredSpliced$readClassType <- res$readClassType
  rowDataFilteredSpliced$novelTranscript <- res$novelTranscript
  end.ptm <- proc.time()
  if (verbose) message("extended annotations for spliced reads in ",
                       round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  return(rowDataFilteredSpliced)
}

.bambu_extract_grl_flat <- function(exonsByReadClass, intronsByReadClass = NULL) {
    ex_cnts <- as.integer(IRanges::elementNROWS(exonsByReadClass))
    unlisted <- unlist(exonsByReadClass, use.names = FALSE)
    ex_mc <- S4Vectors::mcols(unlisted, use.names = FALSE)
    st_chr <- as.character(GenomicRanges::strand(unlisted))
    in_cnts <- integer(length(exonsByReadClass))
    in_s <- integer(0)
    in_e <- integer(0)
    if (!is.null(intronsByReadClass) && length(intronsByReadClass) == length(exonsByReadClass)) {
        in_cnts <- as.integer(IRanges::elementNROWS(intronsByReadClass))
        in_u <- unlist(intronsByReadClass, use.names = FALSE)
        in_s <- as.integer(IRanges::start(in_u))
        in_e <- as.integer(IRanges::end(in_u))
    }
    list(
        seq_names = as.character(Seqinfo::seqlevels(exonsByReadClass)),
        tx_ex_counts = ex_cnts,
        ex_seq_1b = as.integer(GenomicRanges::seqnames(unlisted)),
        ex_start = as.integer(IRanges::start(unlisted)),
        ex_end = as.integer(IRanges::end(unlisted)),
        ex_strand_code = ifelse(st_chr == "+", 1L, ifelse(st_chr == "-", -1L, 0L)),
        ex_rank = as.integer(ex_mc$exon_rank),
        ex_end_rank = as.integer(ex_mc$exon_endRank),
        tx_in_counts = in_cnts,
        in_start = in_s,
        in_end = in_e
    )
}


#' update classificationTable
#' @importFrom GenomicRanges match
#' @noRd
updateWIntronMatches <- function(unlistedIntrons, unlistedIntronsAnnotations,
                                partitioning, classificationTable, annotationGrangesList,
                                rowDataFilteredSpliced, exonsByReadClass, min.exonDistance,
                                min.primarySecondaryDist, min.primarySecondaryDistStartEnd){
  intronMatches <- match(
    unlistedIntrons,unique(unlistedIntronsAnnotations),nomatch = 0) > 0
  intronMatchesList <- relist(intronMatches, partitioning)
  classificationTable$newWithin[all(intronMatchesList) & !(
    classificationTable$compatible == "compatible" |
      classificationTable$equal == "equal")] <- "newWithin"
  intronMatchListRev <- lapply(intronMatchesList, rev)
  lastJunctionMatch <- unlist(lapply(intronMatchListRev,`[[`, 1))
  firstJunctionMatch <- unlist(lapply(intronMatchesList, `[[`, 1))
  strandVec <- rowDataFilteredSpliced$strand
  intronMatchAny <- any(intronMatchesList)
  classificationTable$newLastJunction[which( strandVec == "+" & !lastJunctionMatch &
                                               intronMatchAny | (strandVec == "-" & (!firstJunctionMatch &
                                                                                       intronMatchAny)))] <- "newLastJunction"
  classificationTable$newFirstJunction[which(strandVec == "+" & !firstJunctionMatch &
                                               intronMatchAny | (strandVec == "-" & (!lastJunctionMatch &
                                                                                       intronMatchAny)))] <- "newFirstJunction"
  classificationTable$newJunction[(
    sum(!intronMatchesList) > !firstJunctionMatch + !lastJunctionMatch) &
      intronMatchAny] <- "newJunction"
  classificationTable$allNew[!intronMatchAny] <- "allNew"
  ## assign gene ids based on the max # of matching introns/splice junctions
  overlapsNewIntronsAnnotatedIntrons <- findOverlaps(unlistedIntrons,
                                                     unlistedIntronsAnnotations, type = "equal",
                                                     select = "all", ignore.strand = FALSE)
  if (length(overlapsNewIntronsAnnotatedIntrons)) {
    distNewTxByQuery <- assignGeneIDbyMaxMatch(
      unlistedIntrons,unlistedIntronsAnnotations,
      overlapsNewIntronsAnnotatedIntrons, exonsByReadClass,
      rowDataFilteredSpliced, annotationGrangesList, min.exonDistance,
      min.primarySecondaryDist, min.primarySecondaryDistStartEnd)
    classificationTable$compatible[distNewTxByQuery$queryHits[
      distNewTxByQuery$compatible]] <- "compatible"
    classificationTable$newFirstExon[distNewTxByQuery$queryHits[
      !distNewTxByQuery$startMatch]] <- "newFirstExon"
    classificationTable$newFirstExon[
      classificationTable$newFirstJunction != "newFirstJunction"] <- ""
    classificationTable$newLastExon[distNewTxByQuery$queryHits[
      !distNewTxByQuery$endMatch]] <- "newLastExon"
    classificationTable$newLastExon[
      classificationTable$newLastJunction != "newLastJunction"] <- ""
  }
  return(classificationTable)
}



#' assign gene id by maximum match
#' @importFrom dplyr as_tibble %>% group_by summarise filter ungroup
#' @noRd
assignGeneIDbyMaxMatch <- function(unlistedIntrons,
                                   unlistedIntronsAnnotations, overlapsNewIntronsAnnotatedIntrons,
                                   exonsByReadClass, rowDataFilteredSpliced, annotationGrangesList,
                                   min.exonDistance, min.primarySecondaryDist,
                                   min.primarySecondaryDistStartEnd) {
  maxGeneCountPerNewTx <- as_tibble(data.frame(txId =
                                                 names(unlistedIntrons)[queryHits(overlapsNewIntronsAnnotatedIntrons)],
                                               geneId = mcols(unlistedIntronsAnnotations)$GENEID[
                                                 subjectHits(overlapsNewIntronsAnnotatedIntrons)],
                                               stringsAsFactors = FALSE)) %>%
    group_by(txId, geneId) %>%
    summarise(geneCount = n()) %>%
    group_by(txId) %>%
    filter(geneCount == max(geneCount)) %>%
    filter(!duplicated(txId)) %>%
    ungroup()
  geneIdByIntron <- rep(NA, length(exonsByReadClass))
  geneIdByIntron <- maxGeneCountPerNewTx$geneId[
    match(names(exonsByReadClass), maxGeneCountPerNewTx$txId)]
  rowDataFilteredSpliced$GENEID[is.na(rowDataFilteredSpliced$GENEID)] <-
    geneIdByIntron[is.na(rowDataFilteredSpliced$GENEID)]
  distNewTx <- calculateDistToAnnotation(
    exByTx = exonsByReadClass,
    exByTxRef = annotationGrangesList,
    maxDist = min.exonDistance,
    primarySecondaryDist = min.primarySecondaryDist,
    primarySecondaryDistStartEnd = min.primarySecondaryDistStartEnd,
    ignore.strand = FALSE)
  distNewTxByQuery <- distNewTx %>%
    group_by(queryHits) %>%
    summarise(
      minDist = min(dist),
      startMatch = any(startMatch),
      endMatch = any(endMatch),
      compatible = any(compatible))
  ## note: here is more information that can be used to filter and annotate!
  return(distNewTxByQuery)
}


#' Calculate distance from read class to annotation
#' @param exByTx exByTx
#' @param exByTxRef exByTxRef
#' @param maxDist defaults to 35
#' @param primarySecondaryDist defaults to 5
#' @param ignore.strand defaults to FALSE
#' @importFrom dplyr ungroup %>%
#' @noRd
calculateDistToAnnotation <- function(exByTx, exByTxRef, maxDist = 35,
                                      primarySecondaryDist = 5, primarySecondaryDistStartEnd = 5,
                                      ignore.strand = FALSE, ref_csr_ptr = NULL, ncore = 1L) {
  if (length(exByTx) == 0L || length(exByTxRef) == 0L) {
      return(tibble::tibble(
          queryHits = integer(0), subjectHits = integer(0),
          compatible = logical(0), equal = logical(0), dist = integer(0),
          startMatch = logical(0), endMatch = logical(0),
          uniqueStartLengthQuery = integer(0), uniqueStartLengthSubject = integer(0),
          uniqueEndLengthQuery = integer(0), uniqueEndLengthSubject = integer(0),
          txNumberFiltered = integer(0),
          readClassId = character(0), annotationTxId = character(0), txid = integer(0)
      ))
  }
  if (is.null(ref_csr_ptr)) {
      ref_csr_ptr <- bambu_next_get_csr(exByTxRef, n_threads = ncore)
  }
  ex_cnts <- as.integer(IRanges::elementNROWS(exByTx))
  unlisted <- unlist(exByTx, use.names = FALSE)
  ex_mc <- S4Vectors::mcols(unlisted, use.names = FALSE)
  st_chr <- as.character(GenomicRanges::strand(unlisted))
  res <- .Call(
      "bambu_next_calc_dist_to_annotation_cpp",
      as.character(Seqinfo::seqlevels(exByTx)),
      ex_cnts,
      as.integer(GenomicRanges::seqnames(unlisted)),
      as.integer(IRanges::start(unlisted)),
      as.integer(IRanges::end(unlisted)),
      ifelse(st_chr == "+", 1L, ifelse(st_chr == "-", -1L, 0L)),
      as.integer(ex_mc$exon_rank),
      as.integer(ex_mc$exon_endRank),
      ref_csr_ptr,
      as.integer(maxDist),
      as.integer(primarySecondaryDist),
      as.integer(primarySecondaryDistStartEnd),
      isTRUE(ignore.strand),
      as.integer(ncore),
      PACKAGE = "bambu"
  )
  qh <- res$queryHits
  tx_nf <- if (length(qh) > 0L) as.integer(ave(qh, qh, FUN = length)) else integer(0)
  out <- tibble::tibble(
      queryHits = qh,
      subjectHits = res$subjectHits,
      compatible = res$compatible,
      equal = res$equal,
      dist = res$dist,
      startMatch = res$startMatch,
      endMatch = res$endMatch,
      uniqueStartLengthQuery = res$uniqueStartLengthQuery,
      uniqueStartLengthSubject = res$uniqueStartLengthSubject,
      uniqueEndLengthQuery = res$uniqueEndLengthQuery,
      uniqueEndLengthSubject = res$uniqueEndLengthSubject,
      txNumberFiltered = tx_nf,
      readClassId = names(exByTx)[qh],
      annotationTxId = names(exByTxRef)[res$subjectHits],
      txid = S4Vectors::mcols(exByTxRef)$txid[res$subjectHits]
  )
  return(out)
}


#' generate filtered annotation table
#' @param spliceOverlaps an output from  findSpliceOverlapsByDist()
#' @param primarySecondaryDist default 5
#' @param primarySecondaryDistStartEnd default 5
#' @param exByTx default NULL
#' @param setTMP default NULL
#' @param DistCalculated default FALSE
#' @importFrom dplyr as_tibble group_by %>% mutate n arrange filter arrange
#' @noRd
genFilteredAnTable <- function(spliceOverlaps, primarySecondaryDist = 5,
                               primarySecondaryDistStartEnd = 5, exByTx = NULL, setTMP = NULL,
                               DistCalculated = FALSE) {
  ## initiate the table
  if (isFALSE(DistCalculated)) {
    txToAnTable <- as_tibble(spliceOverlaps) %>% group_by(queryHits) %>%
      mutate(dist = uniqueLengthQuery + uniqueLengthSubject) %>%
      mutate(txNumber = n())
  } else {
    txToAnTable <- as_tibble(spliceOverlaps) %>% group_by(queryHits) %>%
      mutate(dist = uniqueLengthQuery + uniqueLengthSubject +
               uniqueStartLengthQuery + uniqueEndLengthQuery) %>%
      mutate(txNumber = n())
  }
  ## change query hits for step 2 and 3
  if (!is.null(exByTx)) {
    txToAnTable$queryHits <-
      (seq_along(exByTx))[-setTMP][txToAnTable$queryHits]
  }
  ## todo: check filters, what happens to reads with only start and end match?
  if (isFALSE(DistCalculated)) {
    txToAnTableFiltered <- txToAnTable %>%
      group_by(queryHits) %>%
      arrange(queryHits, dist) %>%
      filter(dist <= (min(dist) + primarySecondaryDist)) %>%
      filter(queryElementsOutsideMaxDist + 
               subjectElementsOutsideMaxDist == 
               min(queryElementsOutsideMaxDist +
                     subjectElementsOutsideMaxDist)) %>% 
      filter((uniqueStartLengthQuery <= primarySecondaryDistStartEnd &
                uniqueEndLengthQuery <= primarySecondaryDistStartEnd) ==
               max(uniqueStartLengthQuery <=
                     primarySecondaryDistStartEnd & uniqueEndLengthQuery <=
                     primarySecondaryDistStartEnd)) %>%
      mutate(txNumberFiltered = n())
  } else {
    txToAnTableFiltered <- txToAnTable %>%
      group_by(queryHits) %>%
      arrange(queryHits, dist) %>%
      filter(dist <= (min(dist) + primarySecondaryDist)) %>%
      mutate(txNumberFiltered = n())
  }
  return(txToAnTableFiltered)
}

#' extended annotations for unspliced reads
#' @param se a summarized experient object
#' @param seFilteredSpliced seFilteredSpliced
#' @param exonsByReadClass exonsByReadClass
#' @param annotationGrangesList annotationGrangesList
#' @param filterSet1 filterSet1
#' @param min.exonOverlap min.exonOverlap
#' @param verbose verbose
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom SummarizedExperiment rbind
#' @noRd
addNewUnsplicedReadClasses <- function(rowDataFilteredUnspliced,  rowDataFilteredSpliced,
                                       exonsByReadClass, annotationGrangesList, min.exonOverlap, verbose,
                                       ref_csr_ptr = NULL) {
  start.ptm <- proc.time()
  rowDataCombined <- rowDataFilteredSpliced
  exonRangesCombined <- exonsByReadClass
  names(exonRangesCombined) <- seq_along(exonRangesCombined)
  if (nrow(rowDataFilteredUnspliced)) {
    q_chrs <- as.character(rowDataFilteredUnspliced$chr)
    q_lvls <- if (is.factor(rowDataFilteredUnspliced$chr)) levels(rowDataFilteredUnspliced$chr) else unique(q_chrs)
    q_seq_1b <- match(q_chrs, q_lvls)
    st_chr <- as.character(rowDataFilteredUnspliced$strand)
    q_st_code <- ifelse(st_chr == "+", 1L, ifelse(st_chr == "-", -1L, 0L))

    if (length(annotationGrangesList) > 0L) {
        if (is.null(ref_csr_ptr)) ref_csr_ptr <- bambu_next_get_csr(annotationGrangesList)
        naOverlapUnspliced <- .Call(
            "bambu_next_filter_unspliced_overlap_cpp",
            ref_csr_ptr,
            as.character(q_lvls),
            as.integer(q_seq_1b),
            as.integer(rowDataFilteredUnspliced$start),
            as.integer(rowDataFilteredUnspliced$end),
            as.integer(q_st_code),
            as.integer(min.exonOverlap),
            PACKAGE = "bambu"
        )
    } else {
        naOverlapUnspliced <- seq_len(nrow(rowDataFilteredUnspliced))
    }

    if (length(naOverlapUnspliced)) {
      rowDataFilteredUnspliced <- rowDataFilteredUnspliced[naOverlapUnspliced,]
      rowDataFilteredUnspliced$GENEID <- NA
      rowDataFilteredUnspliced$TXNAME <- NA
      rowDataFilteredUnspliced$readClassType <- "unsplicedNew"
      rowDataFilteredUnspliced$novelGene <- TRUE
      rowDataFilteredUnspliced$novelTranscript <- TRUE

      all_lvls <- unique(c(q_lvls, seqlevels(annotationGrangesList)))
      keep_chrs <- factor(as.character(rowDataFilteredUnspliced$chr), levels = all_lvls)
      exonsByReadClassUnspliced <- GRanges(
        seqnames = keep_chrs,
        ranges = IRanges(start = rowDataFilteredUnspliced$start,
                         end = rowDataFilteredUnspliced$end),
        strand = rowDataFilteredUnspliced$strand)
      partitioning <- PartitioningByEnd(seq_along(exonsByReadClassUnspliced), names = NULL)
      exonsByReadClassUnspliced$exon_rank <- rep(1, length(exonsByReadClassUnspliced))
      exonsByReadClassUnspliced$exon_endRank <- rep(1, length(exonsByReadClassUnspliced))
      exonsByReadClassUnspliced <- relist(exonsByReadClassUnspliced, partitioning)

      rowDataCombined <- bind_rows(rowDataFilteredSpliced, rowDataFilteredUnspliced)
      exonRangesCombined <- c(exonsByReadClass, exonsByReadClassUnspliced)
      names(exonRangesCombined) <- seq_along(exonRangesCombined)
    }
  } 
  end.ptm <- proc.time()
  if (verbose) message("extended annotations for unspliced reads in ",
                       round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  return(list("rowDataCombined" = rowDataCombined, 
              "exonRangesCombined" = exonRangesCombined))
}

#' calculate distance between first and last exon matches
#' @param candidateList candidateList
#' @param filteredOverlapList filteredOverlapList
#' @importFrom dplyr select rename %>% left_join group_by filter 
#'     ungroup distinct
#' @noRd
includeOverlapReadClass <- function(candidateList, filteredOverlapList) {
  temp <- left_join(candidateList, filteredOverlapList,
                    by = c("subjectHits" = "queryHits")) %>%
    group_by(queryHits) %>%
    filter(!subjectHits.y %in% subjectHits, !is.na(subjectHits.y)) %>%
    ungroup() %>%
    dplyr::select(queryHits, subjectHits.y) %>%
    distinct() %>%
    rename(subjectHits = subjectHits.y)
  return(temp)
}

#' combine annotations with predicted transcripts
#' @noRd
combineWithAnnotations <- function(rowDataCombinedFiltered, 
                                        exonRangesCombinedFiltered, annotationGrangesList, prefix,
                                        ref_csr_ptr = NULL, ncore = 1L){
    rowDataCombinedFiltered$txClassDescription <- rowDataCombinedFiltered$readClassType
    rowDataCombinedFiltered$txClassDescription[rowDataCombinedFiltered$readClassType
                                       == "unsplicedNew" & rowDataCombinedFiltered$novelGene] <-
      "newGene-unspliced"
    rowDataCombinedFiltered$txClassDescription[rowDataCombinedFiltered$readClassType
                                       == "allNew" & rowDataCombinedFiltered$novelGene] <-
      "newGene-spliced"
    extendedAnnotationRanges <- exonRangesCombinedFiltered
    mcols(extendedAnnotationRanges) <-
      rowDataCombinedFiltered[, c("GENEID", "novelGene", "novelTranscript", "txClassDescription","readCount", "NDR",
                                  "relReadCount", "relSubsetCount")]
    equalRanges = rowDataCombinedFiltered[!is.na(rowDataCombinedFiltered$TXNAME),]
    #remove extended ranges that are already present in annotation
    extendedAnnotationRanges <- extendedAnnotationRanges[is.na(rowDataCombinedFiltered$TXNAME)]
    annotationRangesToMerge <- annotationGrangesList
    if(length(annotationGrangesList)){
        mcols(annotationRangesToMerge)$readCount <- NA
        mcols(annotationRangesToMerge)$txClassDescription <- "annotation"
        mcols(annotationRangesToMerge)$novelTranscript <- FALSE
        mcols(annotationRangesToMerge)$novelGene <- FALSE
        mcols(annotationRangesToMerge)$NDR <- NA
        mcols(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)[,colnames(mcols(extendedAnnotationRanges))]
        #copy over stats to annotations from read classes
        mcols(annotationRangesToMerge[equalRanges$TXNAME])$NDR = equalRanges$NDR
        mcols(annotationRangesToMerge[equalRanges$TXNAME])$readCount = equalRanges$readCount
        mcols(annotationRangesToMerge[equalRanges$TXNAME])$relReadCount = equalRanges$relReadCount
        mcols(annotationRangesToMerge[equalRanges$TXNAME])$relSubsetCount = equalRanges$relSubsetCount
    }
    n_novel <- length(extendedAnnotationRanges)
    if (n_novel > 0L) {
      mcols(extendedAnnotationRanges)$TXNAME <- paste0(prefix, "Tx", seq_len(n_novel))
      names(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)$TXNAME
      extendedAnnotationRanges <- c(extendedAnnotationRanges, annotationRangesToMerge)
      mcols(extendedAnnotationRanges)$txid <- seq_along(extendedAnnotationRanges)

      if (is.null(ref_csr_ptr) && length(annotationGrangesList) > 0L) {
          ref_csr_ptr <- bambu_next_get_csr(annotationGrangesList, n_threads = ncore)
      }
      ex_counts <- as.integer(IRanges::elementNROWS(extendedAnnotationRanges))
      unlisted <- unlist(extendedAnnotationRanges, use.names = FALSE)
      ex_starts <- as.integer(unlisted@ranges@start)
      ex_ends <- as.integer(ex_starts + unlisted@ranges@width - 1L)
      ex_mc <- S4Vectors::mcols(unlisted, use.names = FALSE)
      ex_rank <- if (!is.null(ex_mc$exon_rank)) as.integer(ex_mc$exon_rank) else seq_along(ex_starts)
      ex_end_rank <- if (!is.null(ex_mc$exon_endRank)) as.integer(ex_mc$exon_endRank) else ex_rank
      first_idx <- cumsum(c(1L, ex_counts[-length(ex_counts)]))
      seq_names <- as.character(Seqinfo::seqlevels(extendedAnnotationRanges))
      seq_lens <- as.integer(Seqinfo::seqlengths(extendedAnnotationRanges)[seq_names])
      seq_lens[is.na(seq_lens)] <- -1L
      tx_seq_1b <- as.integer(GenomicRanges::seqnames(unlisted))[first_idx]
      st_vec <- as.character(GenomicRanges::strand(unlisted))[first_idx]
      tx_strand_code <- ifelse(st_vec == "+", 1L, ifelse(st_vec == "-", -1L, 0L))
      mc <- S4Vectors::mcols(extendedAnnotationRanges, use.names = FALSE)
      tx_names <- as.character(mc$TXNAME)
      tx_gene_names <- as.character(mc$GENEID)

      res <- .Call(
          "bambu_next_extend_min_eq_csr_cpp",
          ref_csr_ptr, as.integer(n_novel),
          seq_names, seq_lens, tx_names, tx_gene_names,
          tx_seq_1b, as.integer(tx_strand_code), ex_counts,
          ex_starts, ex_ends, ex_rank, ex_end_rank, as.integer(ncore),
          PACKAGE = "bambu"
      )
      part <- IRanges::PartitioningByEnd(cumsum(res$eq_lens), names = NULL)
      mcols(extendedAnnotationRanges)$eqClassById <- IRanges::relist(res$eq_vals, part)
      nt <- length(tx_names)
      .Call(
          "bambu_next_csr_cache_put",
          res$csr_ptr, nt, ex_starts, ex_ends, res$eq_lens, res$eq_vals,
          if (nt > 0L) tx_names[1L] else "",
          if (nt > 0L) tx_names[nt] else "",
          PACKAGE = "bambu"
      )
    } else {
      extendedAnnotationRanges <- annotationRangesToMerge
      mcols(extendedAnnotationRanges)$txid <- seq_along(extendedAnnotationRanges)
      mcols(extendedAnnotationRanges)$relReadCount = NA
      mcols(extendedAnnotationRanges)$relSubsetCount = NA
      if (length(extendedAnnotationRanges) > 0L && is.null(mcols(extendedAnnotationRanges)$eqClassById)) {
          extendedAnnotationRanges <- bambu_next_attach_csr(extendedAnnotationRanges, n_threads = ncore, compute_eq = TRUE)
      }
    }
    mcols(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)[, 
                 c("TXNAME", "GENEID", "NDR", "novelGene", "novelTranscript", "txClassDescription","readCount","relReadCount", "relSubsetCount", "txid", "eqClassById")]
    return(extendedAnnotationRanges)
}


#' Estimate distance between read class and annotations
#' @param seReadClass seReadClass
#' @inheritParams bambu
#' @importFrom dplyr select left_join as_tibble mutate %>%
#' @noRd
isore.estimateDistanceToAnnotations <- function(seReadClass,
                                                annotationGrangesList, min.exonDistance = 35,
                                                min.primarySecondaryDist = 5, min.primarySecondaryDistStartEnd = 100000, 
                                                additionalFiltering = FALSE, verbose = FALSE) {
  start.ptm <- proc.time()
  readClassTable <-
    as_tibble(rowData(seReadClass), rownames = "readClassId") %>%
    dplyr::select(readClassId, confidenceType)
  distTable <- calculateDistToAnnotation(rowRanges(seReadClass),
                                         annotationGrangesList, maxDist = min.exonDistance,
                                         primarySecondaryDist = min.primarySecondaryDist,
                                         primarySecondaryDistStartEnd = min.primarySecondaryDistStartEnd,
                                         ignore.strand = FALSE)
  distTable$readCount <- assays(seReadClass)$counts[distTable$readClassId, ] 
  if (additionalFiltering) 
    distTable <- left_join(distTable, select(readClassTable,
                                             readClassId, confidenceType), by = "readClassId") %>%
    mutate(relativeReadCount = readCount / txNumberFiltered)
  distTable <- dplyr::select(distTable, annotationTxId, txid, readClassId,
      readCount, compatible, equal,dist)
  distTable <- left_join(distTable, as_tibble(mcols(annotationGrangesList)[, c("txid", "GENEID")]),
      by = c("txid" = "txid"))
  
  end.ptm <- proc.time()
  if (verbose) message("calculated distance table in ",
                       round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  readClassTable <- addGeneIdsToReadClassTable(readClassTable, distTable, 
                                               seReadClass, verbose)
  distTable <- DataFrame(distTable)
  distTable$eqClassById <- as.list(mcols(annotationGrangesList)$eqClassById)[distTable$txid]
  
  metadata(seReadClass) <- list(distTable = distTable)
  rowData(seReadClass) <- readClassTable
  return(seReadClass)
}



#' generate readClassTable
#' @importFrom dplyr select filter distinct unlist group_by mutate %>% ungroup
#'     left_join
#' @noRd
addGeneIdsToReadClassTable <- function(readClassTable, distTable, 
                                       seReadClass, verbose){
  readClassTable$equal <- readClassTable$readClassId %in%
    unlist((filter(distTable, equal) %>%
              dplyr::select(readClassId) %>% distinct()))
  readClassTable$compatible <- readClassTable$readClassId %in%
    unlist((filter(distTable, compatible) %>% distinct()))
  start.ptm <- proc.time()
  # assign read classes to genes based on the highest read count per gene
  readClassToGeneIdTable <- select(distTable, readClassId, GENEID,
                                   readCount) %>% group_by(GENEID) %>%
    mutate(geneCount = sum(readCount)) %>% distinct() %>%
    group_by(readClassId) %>% filter(geneCount == max(geneCount)) %>%
    filter(row_number() == 1) %>% 
    dplyr::select(readClassId, geneId = GENEID) %>% ungroup()
  newGeneCandidates <- (!readClassTable$readClassId %in%
                          readClassToGeneIdTable$readClassId)
  readClassToGeneIdTableNew <- 
      assignGeneIdsNoReference(rowRanges(seReadClass)[newGeneCandidates], 
      prefix = "unassigned")
  readClassToGeneIdTableNew <- tibble(names(seReadClass)[newGeneCandidates], 
      readClassToGeneIdTableNew)
  colnames(readClassToGeneIdTableNew) <- c('readClassId', 'geneId')
  readClassGeneTable <- rbind(readClassToGeneIdTable, 
                              readClassToGeneIdTableNew)
  readClassTable <- left_join(readClassTable, readClassGeneTable,
                              by = "readClassId") %>%
    dplyr::select(confidenceType, geneId, compatible, equal)
  end.ptm <- proc.time()
  if (verbose) message("added gene Ids for each read class ",
                       round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  return(readClassTable)
}
