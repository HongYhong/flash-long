#' process reads
#' @param reads path to BAM file(s)
#' @param annotations path to GTF file or TxDb object
#' @param genomeSequence path to FA file or BSgenome object
#' @param readClass.outputDir path to readClass output directory
#' @param yieldSize yieldSize
#' @param bpParameters BioParallel parameter
#' @param stranded stranded
#' @param verbose verbose
#' @importFrom Rsamtools yieldSize BamFileList yieldSize<-
#' @importFrom methods is 
#' @importFrom BiocParallel bplapply
#' @importFrom BiocGenerics basename
#' @noRd
bambu.processReads <- function(reads, annotations, genomeSequence,
    readClass.outputDir=NULL, yieldSize=1000000, bpParameters, 
    stranded=FALSE, verbose=FALSE, isoreParameters = setIsoreParameters(NULL),
    lowMemory=FALSE, trackReads = trackReads, fusionMode = fusionMode,
    ncore = NULL) {
    genomeSequence <- checkInputSequence(genomeSequence)
    # ===# create BamFileList object from character #===#
    if (is(reads, "BamFile")) {
        if (!is.null(yieldSize)) {
            yieldSize(reads) <- yieldSize
        } else {
            yieldSize <- yieldSize(reads)
        }
        reads <- BamFileList(reads)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
    } else if (is(reads, "BamFileList")) {
        if (!is.null(yieldSize)) {
            yieldSize(reads) <- yieldSize
        } else {
            yieldSize <- min(yieldSize(reads))
        }
    } else if (any(!grepl("\\.bam$", reads))) {
        stop("Bam file is missing from arguments.")
    } else {
        if (is.null(yieldSize)) yieldSize <- NA
        reads <- BamFileList(reads, yieldSize = yieldSize)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
    }
    min.readCount = isoreParameters[["min.readCount"]]
    fitReadClassModel = isoreParameters[["fitReadClassModel"]]
    defaultModels = isoreParameters[["defaultModels"]]
    returnModel = isoreParameters[["returnModel"]]
    min.exonOverlap = isoreParameters[["min.exonOverlap"]]

    n_threads <- if (!is.null(ncore)) {
        max(1L, as.integer(ncore))
    } else if (!missing(bpParameters) && !is.null(bpParameters$workers)) {
        max(1L, as.integer(bpParameters$workers))
    } else {
        1L
    }
    env_thr <- suppressWarnings(as.integer(Sys.getenv("BAMBU_NEXT_STAGE1_THREADS", "")))
    if (!is.na(env_thr) && env_thr >= 1L) n_threads <- env_thr

    # Pre-build FlatCSRAnnotation once in main process so multi-sample Stage 1 uses zero CoW forks
    if (length(annotations) > 0L) {
        annotations <- bambu_next_attach_csr(annotations, n_threads = n_threads, compute_eq = FALSE)
    }

    readClassList <- lapply(names(reads), function(bamFileName) {
        bambu.processReadsByFile(bam.file = reads[bamFileName],
        genomeSequence = genomeSequence, annotations = annotations,
        readClass.outputDir = readClass.outputDir,
        stranded = stranded, min.readCount = min.readCount,
        fitReadClassModel = fitReadClassModel, min.exonOverlap = min.exonOverlap,
        defaultModels = defaultModels, returnModel = returnModel, verbose = verbose,
        lowMemory = lowMemory, trackReads = trackReads, fusionMode = fusionMode,
        n_threads = n_threads)
    })
    return(readClassList)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @importFrom GenomeInfoDb seqlevels seqlevels<- keepSeqlevels
#' @noRd
bambu.processReadsByFile <- function(bam.file, genomeSequence, annotations,
    readClass.outputDir = NULL, stranded = FALSE, min.readCount = 2, 
    fitReadClassModel = TRUE, min.exonOverlap = 10, defaultModels = NULL, returnModel = FALSE, 
    verbose = FALSE, lowMemory = FALSE, trackReads = FALSE, fusionMode = FALSE,
    n_threads = 1L) {
    if(verbose) message(names(bam.file)[1])
    se <- bambu_next_construct_rc_se_cpp(
        bam_file = bam.file,
        genomeSequence = genomeSequence,
        annotations = annotations,
        stranded = stranded,
        verbose = verbose,
        trackReads = trackReads,
        n_threads = n_threads,
        run_name = names(bam.file)[1L],
        lowMemory = lowMemory
    )
    se <- scoreReadClasses(se, genomeSequence, annotations,
                             defaultModels = defaultModels,
                             fit = fitReadClassModel,
                             returnModel = returnModel,
                             min.readCount = min.readCount,
                             min.exonOverlap = min.exonOverlap,
                             fusionMode = fusionMode,
                             verbose = verbose,
                             n_threads = n_threads)
    if (!is.null(readClass.outputDir)) {
        readClassFile <- paste0(readClass.outputDir,names(bam.file),
                                "_readClassSe.rds")
        if (file.exists(readClassFile)) {
            show(paste(readClassFile, "exists, will be overwritten"))
            warning(readClassFile, "exists, will be overwritten")
        } else {
            readClassFile <- BiocFileCache::bfcnew(BiocFileCache::BiocFileCache(
                readClass.outputDir, ask = FALSE),
                paste0(names(bam.file),"_readClassSe"), ext = ".rds")
        }
        saveRDS(se, file = readClassFile)
        se <- readClassFile
    }
    return(se)
}

lowMemoryConstructReadClasses <- function(readGrgList, genomeSequence, 
                                          annotations, stranded, verbose,bam.file){
    readGrgList = split(readGrgList, getChrFromGrList(readGrgList))
    se = lapply(names(readGrgList),FUN = function(i){
        if(length(readGrgList[[i]]) == 0) return(NULL)
        # create error and strand corrected junction tables
        unlisted_junctions <- unlistIntrons(readGrgList[[i]], use.ids = TRUE)
        uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
                                                         annotations,genomeSequence, stranded = stranded, verbose = verbose)
        se.temp <- isore.constructReadClasses(readGrgList[[i]], 
                                              unlisted_junctions, uniqueJunctions, runName = names(bam.file)[1],
                                              annotations, stranded, verbose)
        return(se.temp)
    })
    se = se[!sapply(se, FUN = is.null)]
    se = do.call("rbind",se)
    rownames(se) = paste("rc", seq_len(nrow(se)), sep = ".")
    return(se)
}

#' Check seqlevels for reads and annotations
#' @importFrom GenomeInfoDb seqlevels
#' @noRd
seqlevelCheckReadsAnnotation <- function(reads, annotations){
    warnings = c()
    if (length(intersect(seqlevels(reads),
                         seqlevels(annotations))) == 0)
        warnings = c(warnings, paste0("no annotations with matching seqlevel styles, ",
        "all missing chromosomes will use de-novo annotations"))
    if (!all(seqlevels(reads) %in% 
             seqlevels(annotations))) 
        warnings = c(warnings, paste0("not all chromosomes present in reference annotations, ",
            "annotations might be incomplete. Please compare objects ",
            "on the same reference"))
    return(warnings)
}