#' Reduce transcript expression to gene expression
#' @title transcript to gene expression
#' @param se a summarizedExperiment object from \code{\link{bambu}}
#' @return A SummarizedExperiment object
#' @import data.table 
#' @export
#' @examples
#' se <- readRDS(system.file("extdata",
#'     "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' transcriptToGeneExpression(se)
transcriptToGeneExpression <- function(se) {
    cnt_mat <- as.matrix(assays(se)$counts)
    runnames <- rename_duplicatedNames(colnames(cnt_mat))
    colnames(cnt_mat) <- runnames
    colData(se)@rownames <- rename_duplicatedNames(colData(se)@rownames)

    rowDataSe <- as.data.table(rowData(se))
    tx_genes <- as.character(rowDataSe$GENEID)
    tx_gene_sum <- rowsum(cnt_mat, group = tx_genes, reorder = FALSE, na.rm = TRUE)

    incomp <- data.table::copy(metadata(se)$incompatibleCounts)
    if (!is.null(incomp) && nrow(incomp) > 0L) {
        incomp_cols <- intersect(runnames, setdiff(colnames(incomp), c("GENEID", "TXNAME")))
        all_genes <- sort(unique(c(rownames(tx_gene_sum), as.character(incomp$GENEID))), method = "radix")
        gene_mat <- matrix(0.0, nrow = length(all_genes), ncol = length(runnames),
                           dimnames = list(all_genes, runnames))
        m_tx <- match(rownames(tx_gene_sum), all_genes)
        gene_mat[m_tx, ] <- tx_gene_sum
        if (length(incomp_cols) > 0L) {
            incomp_vals <- as.matrix(incomp[, ..incomp_cols])
            if (anyDuplicated(incomp$GENEID)) {
                incomp_vals <- rowsum(incomp_vals, group = as.character(incomp$GENEID), reorder = FALSE, na.rm = TRUE)
                incomp_genes <- rownames(incomp_vals)
            } else {
                incomp_genes <- as.character(incomp$GENEID)
            }
            m_inc <- match(incomp_genes, all_genes)
            c_inc <- match(incomp_cols, runnames)
            gene_mat[m_inc, c_inc] <- gene_mat[m_inc, c_inc] + incomp_vals
        }
    } else {
        all_genes <- sort(rownames(tx_gene_sum), method = "radix")
        gene_mat <- tx_gene_sum[all_genes, , drop = FALSE]
    }

    col_totals <- pmax(colSums(gene_mat), 1.0)
    cpm_mat <- sweep(gene_mat, 2L, col_totals, "/") * 1e6

    ## geneRanges
    exByGene <- reducedRangesByGenes(rowRanges(se))
    if ("txClassDescription" %in% colnames(rowDataSe)) {
        rowDataSe <- rowDataSe[, .(TXNAME, GENEID, txClassDescription)]
        rowDataSe[, newGeneClass := ifelse(grepl("ENSG", GENEID),
            "annotation", unique(txClassDescription)), by = GENEID]
        mcols(exByGene) <- unique(rowDataSe[, .(GENEID,
            newGeneClass)])[match(names(exByGene), GENEID)]
    }

    ColData <- colData(se)
    ColData@rownames <- runnames
    ColData@listData$name <- runnames
    seOutput <- SummarizedExperiment(
        assays = SimpleList(counts = gene_mat, CPM = cpm_mat),
        rowRanges = exByGene[all_genes],
        colData = ColData
    )
    return(seOutput)
}