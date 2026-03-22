#' @title Publication & disease–gene query tools
#' @description Search PubMed and DisGeNET for term overlaps and counts.
#' @name publication_analysis
NULL

#' Perform a PubMed AND query of two terms
#'
#' @param query1,query2 Character; individual PubMed queries
#' @return A PubMed ID object from easyPubMed::get_pubmed_ids()
#' @import easyPubMed
#' @export
findDoubleQueryPapers <- function(query1, query2) {
  combined <- paste0("(", query1, ") AND (", query2, ")")
  easyPubMed::get_pubmed_ids(combined)
}

#' Get counts for a term and term+autophagy overlap
#'
#' @param query Character; PubMed query
#' @param statsOnly Logical; if TRUE, return counts only
#' @return Numeric vector: c(totalCount, autophagyCount) or full PubMed objects
#' @export
findAutophagyDoublePapers <- function(query, statsOnly = TRUE) {
  ids1 <- easyPubMed::get_pubmed_ids(query)
  ids2 <- findDoubleQueryPapers(query, "autophagy")
  if (statsOnly) {
    return(c(as.numeric(ids1$Count), as.numeric(ids2$Count)))
  }
  message("Total: ", ids1$Count, "; Autophagy: ", ids2$Count)
  list(all = ids1, autophagy = ids2)
}

#' Build a table of autophagy stats for multiple terms
#'
#' @param term Character; a single query term
#' @param table List of terms
#' @return Matrix of counts for each term
#' @export
getTermPapers <- function(term, table) {
  t(sapply(table, findAutophagyDoublePapers, statsOnly = TRUE))
}

#' Count autophagy papers per gene
#'
#' @param input Character; gene ID
#' @param gene2pubmed_data Two-column matrix/data.frame gene→PubMed ID
#' @param autophagy_data Character vector of autophagy PubMed IDs
#' @return Integer vector c(totalPapers, autophagyPapers)
#' @export
autophagyPapersOneGene <- function(input, gene2pubmed_data, autophagy_data) {
  papers <- gene2pubmed_data[gene2pubmed_data[,1] == input, 2]
  c(length(papers), length(intersect(papers, autophagy_data)))
}

#' Apply autophagyPapersOneGene across all unique genes
#'
#' @param gene2pubmed_data Two-column matrix/data.frame gene→PubMed ID
#' @param autophagy_data Character vector of autophagy PubMed IDs
#' @return List of c(total, autophagy) per gene
#' @export
autophagyPapersAllGenes <- function(gene2pubmed_data, autophagy_data) {
  genes <- unique(gene2pubmed_data[,1])
  sapply(genes, autophagyPapersOneGene,
         gene2pubmed_data = gene2pubmed_data,
         autophagy_data    = autophagy_data)
}
