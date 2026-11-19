#' @title Write Bambu results to GTF and transcript/gene-count files
#' @param se a \code{\link{SummarizedExperiment}} object 
#' from \code{\link{bambu}}.
#' @param path the destination of the output files 
#' (gtf, transcript counts, and gene counts)
#' @param prefix the prefix of the output files
#' @details The function will write the output from Bambu to files. The 
#' annotations will be written to a .gtf file, transcript counts (total counts, 
#' CPM, full-length counts, partial-length counts, and unique counts) and gene counts
#' will be written to .txt files.  
#' @export
#' @examples
#' se <- readRDS(system.file("extdata",
#'     "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' path <- tempdir()
#' writeBambuOutput(se, path)
writeBambuOutput <- function(se, path, prefix = "") {
    if (missing(se) | missing(path)) {
        stop("Both summarizedExperiment object from bambu and
            the path for the output files are required.")
    } else {
        outdir <- paste0(path, "/")
        if (!dir.exists(outdir))
            dir.create(outdir, recursive = TRUE)

        transcript_grList <- rowRanges(se)
        transcript_gtffn <- paste(outdir, prefix,
            "extended_annotations.gtf", sep = "")
        gtf <- writeToGTF(annotation = transcript_grList,
            file = transcript_gtffn)
        
        for(d in names(assays(se))){
            writeCountsOutput(se, varname=d,
                             feature='transcript',outdir, prefix)
        }
        seGene <- transcriptToGeneExpression(se)
        writeCountsOutput(seGene, varname='counts',
                             feature='gene',outdir, prefix)
    }
}


#' helper function to write counts 
#' @noRd
writeCountsOutput <- function(se, varname = "counts",
                              feature = "transcript", outdir, prefix){
    estimates <- data.table(as.data.frame(assays(se)[[varname]]),
                            keep.rownames = TRUE) 
    if(feature == "transcript"){
        setnames(estimates, "rn", "TXNAME")
        geneIDs <- data.table(as.data.frame(rowData(se))[,c("TXNAME","GENEID")])
        estimates <- geneIDs[estimates, on = "TXNAME"]
    }else{
        setnames(estimates, "rn","GENEID")
    }
    estimatesfn <- paste(outdir, prefix, 
                                 varname,"_",feature,".txt", sep = "")
    data.table::fwrite(estimates, file = estimatesfn,
                       sep = "\t", quote = FALSE, row.names = FALSE)
}

#' Write annotation GRangesList into a GTF file
#' @title write GRangeslist into GTF file
#' @param annotation a \code{GRangesList} object
#' @param file the output gtf file name
#' @param geneIDs an optional dataframe of geneIDs (column 2) with
#'     the corresponding transcriptIDs (column 1)
#' @return gtf a GTF dataframe
#' @importFrom dplyr select as_tibble mutate %>% left_join arrange group_by
#'     ungroup recode_factor
#' @importFrom methods is
#' @export
#' @examples
#' outputGtfFile <- tempfile()
#' gr <- readRDS(system.file("extdata",
#'     "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' writeToGTF(gr, outputGtfFile)
writeToGTF <- function(annotation, file, geneIDs = NULL) {
    if (missing(annotation) | missing(file)) {
        stop("Both GRangesList and the name of the output file are required.")
    } else if (!is(annotation, "CompressedGRangesList")) {
        stop("The inputted GRangesList is of the wrong class.")
    }
    tx_names <- as.character(names(annotation))
    if (missing(geneIDs) || is.null(geneIDs)) {
        mc <- S4Vectors::mcols(annotation, use.names = FALSE)
        if (!is.null(mc$GENEID)) {
            if (!is.null(mc$TXNAME)) {
                tx_gene_names <- as.character(mc$GENEID[match(tx_names, as.character(mc$TXNAME))])
            } else {
                tx_gene_names <- as.character(mc$GENEID)
            }
        } else {
            tx_gene_names <- rep("", length(annotation))
        }
    } else {
        tx_gene_names <- as.character(geneIDs$GENEID[match(tx_names, as.character(geneIDs$TXNAME))])
    }

    tx_ex_counts <- as.integer(IRanges::elementNROWS(annotation))
    unlisted <- unlist(annotation, use.names = FALSE)
    ex_mc <- S4Vectors::mcols(unlisted, use.names = FALSE)
    st_chr <- as.character(GenomicRanges::strand(unlisted))
    ex_st_code <- ifelse(st_chr == "+", 1L, ifelse(st_chr == "-", -1L, 0L))
    ex_rank <- if (!is.null(ex_mc$exon_rank)) as.integer(ex_mc$exon_rank) else seq_len(length(unlisted))

    .Call(
        "bambu_next_write_gtf_cpp",
        as.character(path.expand(file)),
        as.character(Seqinfo::seqlevels(annotation)),
        tx_names,
        tx_gene_names,
        tx_ex_counts,
        as.integer(GenomicRanges::seqnames(unlisted)),
        as.integer(IRanges::start(unlisted)),
        as.integer(IRanges::end(unlisted)),
        as.integer(ex_st_code),
        as.integer(ex_rank),
        PACKAGE = "bambu"
    )
    invisible(NULL)
}


#' Outputs GRangesList object from reading a GTF file
#' @title convert a GTF file into a GRangesList
#' @param file a .gtf file
#' @param keep.extra.columns a vector with names of columns to keep from 
#' the the attributes in the gtf file. For ensembl, this could be 
#' keep.extra.columns=c('gene_name','gene_biotype',
#' 'transcript_biotype', 'transcript_name')
#' @return grlist a \code{GRangesList} object, with two columns
#' \describe{
#'     \item{TXNAME}{specifying prefix for new gene Ids (genePrefix.number),
#'         defaults to empty}
#'     \item{GENEID}{indicating whether filter to remove read classes which are
#'         a subset of known transcripts(), defaults to TRUE}
#'   }
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @export
#' @examples
#' gtf.file <- system.file("extdata",
#'     "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
#'     package = "bambu"
#' )
#' readFromGTF(gtf.file)
readFromGTF <- function(file, keep.extra.columns = NULL){
    if (missing(file)) {
        stop('A GTF file is required.')
    }else{
        data <- utils::read.delim(file,header = FALSE, comment.char = '#')
        colnames(data) <- c("seqname","source","type","start",
            "end","score","strand","frame","attribute")
    data <- data[data$type == 'exon',]
    data$strand[data$strand == '.'] <- '*'
    data$GENEID = gsub('gene_id (.*?);.*','\\1',data$attribute)
    data$TXNAME = gsub('.*transcript_id (.*?);.*', '\\1',data$attribute)
    data$exon_rank = gsub('.*exon_number (.*?);.*', '\\1',data$attribute)
    if (!is.null(keep.extra.columns)) {
        for (extraColumn in seq_along(keep.extra.columns)) {
            data[,keep.extra.columns[extraColumn]] <-
                gsub(paste0('.*',keep.extra.columns[extraColumn],
                ' (.*?);.*'), '\\1',data$attribute)
            data[grepl(';', data[,keep.extra.columns[extraColumn]]),
                keep.extra.columns[extraColumn]] <- ''
        }
    }
    grlist <- makeGRangesListFromDataFrame(
        data[,c('seqname', 'start','end','strand','TXNAME','exon_rank')],
        split.field = 'TXNAME',keep.extra.columns = TRUE)
    grlist <- grlist[IRanges::order(start(grlist))]
    geneData = (unique(data[,c('TXNAME', 'GENEID',keep.extra.columns)]))
    mcols(grlist) <-
        DataFrame(geneData[(match(names(grlist), geneData$TXNAME)),])
    }
    return(grlist)
}
