#' @title Ensembl BioMart helper functions
#' @description Download gene annotations via biomaRt for human and mouse.
#' @name biomart_helpers
NULL

#' Download Genes Annotated to a GO Term from Ensembl (Human)
#'
#' @param request Character; GO term ID
#' @param text    Logical; if TRUE, write CSV to working dir
#' @param children Logical; if TRUE, include child terms
#' @param attributes Character vector; Ensembl attributes
#' @return data.frame of annotated genes
#' @import biomaRt
#' @export
DownloadGenesInGO <- function(request,
                              text      = FALSE,
                              children  = TRUE,
                              attributes= c("ensembl_gene_id","entrezgene_id","external_gene_name")) {
  ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
  filter <- if (children) "go_parent_term" else "go"
  df <- biomaRt::getBM(attributes = attributes,
                       filters    = filter,
                       values     = request,
                       mart       = ensembl)
  if (text) write.csv(df, "genes_in_go.csv", row.names = FALSE)
  df
}

#' Map a list of identifiers to gene metadata
#'
#' @param input Vector of identifiers.
#' @param filters BioMart filter name.
#' @param attributes Character vector of BioMart attributes to return.
#' @return A data frame of mapped identifiers.
#' @export
getGenesFromList <- function(input = NULL,
                             filters="ensembl_gene_id",
                             attributes= c("ensembl_gene_id","entrezgene_id","external_gene_name")) {
  ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
  biomaRt::getBM(attributes = attributes,
                 filters    = filters,
                 values     = input,
                 mart       = ensembl)
}

#' Download mouse genes annotated to a GO term
#'
#' @param request GO term identifier.
#' @param text Logical; if TRUE, write output to CSV.
#' @param children Logical; if TRUE, query child terms.
#' @param attributes Character vector of BioMart attributes to return.
#' @return A data frame of genes.
#' @export
DownloadGenesInMouseGO <- function(request,
                                   text      = FALSE,
                                   children  = TRUE,
                                   attributes= c("ensembl_gene_id","entrezgene_id","external_gene_name")) {
  ensembl <- biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl")
  filter <- if (children) "go_parent_term" else "go"
  df <- biomaRt::getBM(attributes = attributes,
                       filters    = filter,
                       values     = request,
                       mart       = ensembl)
  if (text) write.csv(df, "mouse_genes_in_go.csv", row.names = FALSE)
  df
}
