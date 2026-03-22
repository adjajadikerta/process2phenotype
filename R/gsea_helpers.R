#' Run GSEA against DisGeNET gene sets
#'
#' @param scores Named numeric vector of gene-level scores.
#' @return A GSEA result object.
#' @importFrom DOSE gseDGN
#' @importFrom clusterProfiler setReadable
#' @import org.Hs.eg.db
#' @export
getGSEA_DGN <- function(scores) {
  start <- Sys.time()
  geneList <- sort(scores, decreasing = TRUE)
  res <- DOSE::gseDGN(
    geneList,
    minGSSize = 80,
    maxGSSize = 90,
    pvalueCutoff = 1,
    pAdjustMethod = "BH"
  )
  res <- clusterProfiler::setReadable(res, "org.Hs.eg.db")
  message("Elapsed: ", Sys.time() - start)
  res
}

#' Run GSEA against Disease Ontology gene sets
#'
#' @param scores Named numeric vector of gene-level scores.
#' @return A GSEA result object.
#' @importFrom clusterProfiler gseDO
#' @export
getGSEA_DO <- function(scores) {
  geneList <- sort(scores, decreasing = TRUE)
  clusterProfiler::gseDO(
    geneList,
    minGSSize = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = TRUE,
    eps = 0
  )
}

#' Run GSEA against MeSH gene sets
#'
#' @param scores Named numeric vector of gene-level scores.
#' @return A GSEA result object.
#' @importFrom meshes gseMeSH
#' @export
getGSEA_MeSH <- function(scores) {
  geneList <- sort(scores, decreasing = TRUE)
  meshes::gseMeSH(
    geneList,
    MeSHDb = "MeSH.Hsa.eg.db",
    database = "gene2pubmed",
    category = "C",
    eps = 0,
    pvalueCutoff = 1
  )
}

#' Run GSEA on a custom TERM2GENE table
#'
#' @param scores Named numeric vector of gene-level scores.
#' @param term2gene A data frame with TERM and GENE columns.
#' @param term2name Optional TERM-to-name mapping.
#' @return A GSEA result object.
#' @importFrom clusterProfiler GSEA
#' @export
getGSEA_custom <- function(scores, term2gene, term2name = NA) {
  start <- Sys.time()
  geneList <- sort(scores, decreasing = TRUE)
  res <- clusterProfiler::GSEA(
    geneList,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 1,
    minGSSize = 5,
    pAdjustMethod = "BH",
    eps = 0
  )
  message("Elapsed: ", Sys.time() - start)
  res
}