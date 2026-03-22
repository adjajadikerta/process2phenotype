#' @title Genomic locus gene lookup
#' @description Retrieve all genes overlapping specified genomic ranges.
#' @name genomic_helpers
NULL

#' Download Genes Within Genomic Loci via Ensembl
#'
#' @param loci A data.frame with columns chromosome_name, start, end
#' @param text Logical; if TRUE, write TSV to \code{filepath}
#' @param filepath Character; directory or prefix for output file
#' @return data.table of genes in loci
#' @import biomaRt
#' @import data.table
#' @export
downloadGenesInLoci <- function(loci = NULL, text = FALSE, filepath = "") {
  if (is.null(loci) || !all(c("chromosome_name","start","end") %in% colnames(loci))) {
    stop("`loci` must be a data.frame with chromosome_name, start, end")
  }
  ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
  getGenesInRange <- function(range) {
    biomaRt::getBM(
      attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name"),
      filters    = c("chromosome_name","start","end"),
      values     = range,
      mart       = ensembl
    )
  }
  loci_list <- split(as.data.frame(loci), seq_len(nrow(loci)))
  results   <- lapply(loci_list, getGenesInRange)
  combined  <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  if (text) {
    out <- file.path(filepath, "genes_in_loci.tsv")
    data.table::fwrite(combined, out, sep = "\t")
  }
  combined
}
