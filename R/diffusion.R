#' @title Network diffusion and analysis functions
#' @description Tools for diffusing scores over igraph networks, plotting results,
#'   and performing degree-stratified statistical tests.
#' @name diffusion
NULL

#' Diffuse a network from a list of genes
#'
#' @param network igraph object to diffuse over
#' @param inputGenes Character vector of seed gene IDs
#' @param targetGenes Character vector of genes to score (defaults to all vertices)
#' @param met Character; diffusion method (e.g. "raw")
#' @param kernel Numeric; diffusion kernel parameter
#' @return Named numeric vector of diffusion scores
#' @import diffuStats
#' @import igraph
#' @export
diffuseList <- function(network, inputGenes, targetGenes = NULL, met = "raw", kernel) {
  if (is.null(targetGenes)) {
    targetGenes <- igraph::V(network)$name
  }
  if (!all(inputGenes %in% targetGenes)) {
    stop("Input genes not contained within target genes")
  }
  testScores <- setNames((targetGenes %in% inputGenes) * 1, targetGenes)
  diffuStats::diffuse(network, scores = testScores, method = met, K = kernel)
}

#' Two-layer diffusion with optional high-score seeds
#'
#' @inheritParams diffuseList
#' @param inputGenesHighScore Character vector of seeds to upweight
#' @return Named numeric vector of diffusion scores
#' @import diffuStats
#' @import igraph
#' @export
diffuseList2Layer <- function(network, inputGenes, inputGenesHighScore = NULL,
                              targetGenes = NULL, met = "raw", kernel) {
  if (is.null(targetGenes)) {
    targetGenes <- igraph::V(network)$name
  }
  if (!all(inputGenes %in% targetGenes)) {
    stop("Input genes not contained within target genes")
  }
  testScores <- setNames((targetGenes %in% inputGenes) * 1, targetGenes)
  if (!is.null(inputGenesHighScore)) {
    testScores[names(testScores) %in% inputGenesHighScore] <- 5
  }
  diffuStats::diffuse(network, scores = testScores, method = met, K = kernel)
}

#' Boxplot of diffusion scores for target vs non-target
#'
#' @param diffusionScores Named numeric vector of diffusion scores
#' @param network igraph object
#' @param positives Character vector of “positive” gene IDs
#' @param target Character vector of genes to include (defaults to all)
#' @param method Character; diffusion method label
#' @param positiveCategory Character; label for positives
#' @param category2 Character vector of a second category (optional)
#' @param category2name Character; label for second category
#' @param returnData Logical; if TRUE, return the data.frame instead of plotting
#' @return ggplot2 object or data.frame if returnData = TRUE
#' @import ggplot2
#' @import igraph
#' @export
targetBoxPlot <- function(diffusionScores, network, positives,
                          target = NULL, method = NULL,
                          positiveCategory = NULL,
                          category2 = NULL, category2name = "Category 2",
                          returnData = FALSE) {
  if (is.null(method)) method <- "raw"
  if (is.null(target)) {
    target <- igraph::V(network)$name %in% names(diffusionScores)
  }
  if (is.null(positiveCategory)) positiveCategory <- "positive category"
  
  df <- data.frame(
    Protein        = igraph::V(network)$name,
    Class          = ifelse(igraph::V(network)$name %in% positives, positiveCategory, "Other"),
    DiffusionScore = diffusionScores,
    Target         = target,
    Method         = method,
    stringsAsFactors = FALSE
  )
  if (!is.null(category2)) {
    df$Class[df$Protein %in% setdiff(category2, positives)] <- category2name
  }
  if (returnData) {
    return(df)
  }
  ggplot2::ggplot(subset(df, Target), aes(x = Class, y = DiffusionScore)) +
    ggplot2::scale_y_log10(limits = stats::quantile(df$DiffusionScore, c(0.02, 0.995))) +
    ggplot2::xlab("Protein class") +
    ggplot2::ylab("Diffusion score") +
    ggplot2::ggtitle(paste("Target proteins in", positiveCategory))
}

#' Plot ROC curve from scores
#'
#' @param scores Numeric vector of scores
#' @param correct Logical vector of true/false labels
#' @return Invisible; side-effect plot plus printed AUC
#' @import ROCR
#' @export
plotROC <- function(scores, correct) {
  pred <- ROCR::prediction(scores, correct)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  plot(perf); abline(a = 0, b = 1)
  auc  <- ROCR::performance(pred, "auc")@y.values[[1]]
  message("AUC = ", auc)
}

#' Plot Precision–Recall curve
#'
#' @inheritParams plotROC
#' @return Invisible; side-effect plot plus printed AUC (PR)
#' @import ROCR
#' @export
plotPRC <- function(scores, correct) {
  pred   <- ROCR::prediction(scores, correct)
  perf   <- ROCR::performance(pred, "prec", "rec")
  plot(perf)
  aucpr  <- ROCR::performance(pred, "aucpr")@y.values[[1]]
  message("AUC (PR) = ", aucpr)
}

#' Half-known test for diffusion performance
#'
#' Randomly hides half of positives and negatives, diffuses from the rest,
#' then plots boxplot & ROC and returns summary.
#'
#' @param perc Fraction of nodes to use as known (0–1)
#' @param network igraph object
#' @param positives Character vector of positive gene IDs
#' @param method Character; diffusion method
#' @param kernel Numeric; diffusion kernel parameter
#' @param categoryName Character; label for positives in plots
#' @return List with diffusion_scores, node summaries, known summaries, target space, positive_scores
#' @import diffuStats
#' @import igraph
#' @import ggplot2
#' @import ROCR
#' @export
halfKnownTest <- function(perc = .5, network, positives,
                          method = "raw", kernel,
                          categoryName = "category") {
  nodes_A    <- igraph::V(network)$name[igraph::V(network)$name %in% positives]
  nodes_notA <- setdiff(igraph::V(network)$name, nodes_A)
  
  known_A   <- sample(nodes_A,    perc * length(nodes_A))
  known_notA<- sample(nodes_notA, perc * length(nodes_notA))
  known     <- c(known_A, known_notA)
  
  target_A  <- setdiff(nodes_A,    known_A)
  target_notA <- setdiff(nodes_notA, known_notA)
  target    <- c(target_A, target_notA)
  
  scores_positive <- igraph::V(network)$name %in% positives
  scores_A        <- setNames((known %in% known_A)*1, known)
  
  diff <- diffuStats::diffuse(network, scores = scores_A, method = method, K = kernel)
  
  df_plot <- data.frame(Scores = diff, Positive = scores_positive)
  df_target <- subset(df_plot, rownames(df_plot) %in% target)
  
  targetBoxPlot(diffusionScores = diff, network = network,
                positives = nodes_A, target = target,
                method = method, positiveCategory = categoryName)
  plotROC(df_target$Scores, df_target$Positive)
  
  list(
    Diffusion_scores     = diff,
    Nodes_summary        = list(Nodes_in_category = nodes_A, Nodes_not_in_category = nodes_notA),
    Known_input_summary  = list(Known_in_category = known_A, Known_not_in_category = known_notA),
    Target_space         = target,
    Positive_scores      = scores_positive
  )
}

#' Plot the top-scoring subnetwork
#'
#' @param diffusionScores Named numeric vector of diffusion scores
#' @param known_positive Character vector of known positives
#' @param number Integer; how many top nodes to include
#' @param network igraph object
#' @return Invisible; side-effect plot
#' @import igraph
#' @import diffuStats
#' @export
topScoresSubnetwork <- function(diffusionScores, known_positive, number, network) {
  ids    <- head(order(diffusionScores, decreasing = TRUE), number)
  names  <- names(diffusionScores)[ids]
  sub    <- igraph::induced.subgraph(network, vids = ids)
  igraph::plot.igraph(
    sub,
    vertex.color = diffuStats::scores2colours(diffusionScores[names]),
    vertex.shape = diffuStats::scores2shapes(names %in% known_positive),
    vertex.label.color = "gray10",
    main = paste0("Top ", number, " proteins from diffusion scores")
  )
}

#' Discretise network by log-degree
#'
#' @param network igraph object
#' @param method Character; discretization method (e.g. "sturges")
#' @return Data.frame with Gene.ID, Log.Degree, Category
#' @import igraph
#' @import varrank
#' @export
networkDiscretisation <- function(network, method = "sturges") {
  deg   <- igraph::degree(network)
  logdeg<- log(deg)
  cats  <- varrank::discretization(logdeg, discretization.method = method)
  data.frame(Gene.ID    = names(logdeg),
             Log.Degree = logdeg,
             Category   = cats,
             row.names  = NULL,
             stringsAsFactors = FALSE)
}

#' Discretise a numeric variable after log-transform
#'
#' @param variable Numeric vector
#' @param names Character vector of IDs (same length as variable)
#' @param method Character; discretization method
#' @return Data.frame with ID, Variable, LogVariable, VarDiscret, LogVarDiscret
#' @import varrank
#' @export
logVariableDiscretisation <- function(variable, names, method = "sturges") {
  df <- data.frame(ID = names, Variable = variable, stringsAsFactors = FALSE)
  df$Variable[df$Variable == 0] <- NA
  df$Variable[is.na(df$Variable)] <- min(df$Variable, na.rm = TRUE) / 2
  df$LogVariable <- log10(df$Variable)
  cats <- varrank::discretization(df$LogVariable, discretization.method = method)
  data.frame(df, VarDiscret = cats, LogVarDiscret = cats, stringsAsFactors = FALSE)
}

#' Median-based degree-stratified Monte Carlo test
#'
#' @param network igraph object
#' @param sample Character vector of gene IDs to test
#' @param method Character; discretization method
#' @param diffusion.scores Named numeric vector of diffusion scores
#' @param permutations Integer; number of replicates
#' @param plottitle Character; title for density plots
#' @return List with Permutations.Median, Test.Median, Test.Percentile
#' @import igraph
#' @import varrank
#' @export
degreeStratifiedMonteCarlo <- function(network, sample,
                                       method = "sturges",
                                       diffusion.scores,
                                       permutations = 1000,
                                       plottitle = "") {
  net_df <- networkDiscretisation(network, method)
  samp_df <- merge(data.frame(Gene.ID = sample, stringsAsFactors = FALSE),
                   net_df, by = "Gene.ID")
  w1 <- table(samp_df$Category)
  w2 <- table(net_df$Category)
  frac <- as.numeric(w1 / w2[names(w1)])
  names(frac) <- names(w1)
  probs <- frac[net_df$Category]
  mc <- replicate(permutations,
                  median(diffusion.scores[sample(net_df$Gene.ID,
                                                 size = length(sample),
                                                 prob = probs)]))
  test_med <- median(diffusion.scores[sample])
  pct <- stats::ecdf(mc)(test_med)
  plot(stats::density(mc), main = plottitle, xlab = "Diffusion score")
  abline(v = test_med, col = "blue")
  list(Permutations.Median = mc,
       Test.Median        = test_med,
       Test.Percentile    = pct)
}

#' Mean-based degree-stratified Monte Carlo test
#'
#' @inheritParams degreeStratifiedMonteCarlo
#' @return List with Permutations.Mean, Test.Mean, Test.Percentile
#' @import igraph
#' @import varrank
#' @export
degreeStratifiedMonteCarloMean <- function(network, sample,
                                           method = "sturges",
                                           diffusion.scores,
                                           permutations = 1000) {
  net_df <- networkDiscretisation(network, method)
  samp_df<- merge(data.frame(Gene.ID = sample, stringsAsFactors = FALSE),
                  net_df, by = "Gene.ID")
  w1 <- table(samp_df$Category)
  w2 <- table(net_df$Category)
  frac <- as.numeric(w1 / w2[names(w1)])
  names(frac) <- names(w1)
  probs <- frac[net_df$Category]
  mc <- replicate(permutations,
                  mean(diffusion.scores[sample(net_df$Gene.ID,
                                               size = length(sample),
                                               prob = probs)]))
  test_mean <- mean(diffusion.scores[sample])
  pct <- stats::ecdf(mc)(test_mean)
  list(Permutations.Mean    = mc,
       Test.Mean          = test_mean,
       Test.Percentile    = pct)
}

#' Retrieve top-scoring genes with metadata
#'
#' @param number Integer; number of top entries
#' @param list Named numeric vector of scores
#' @return Data.frame of gene metadata plus diffusion_score
#' @export
getTopList <- function(number = 200, list) {
  top_ids <- head(order(list, decreasing = TRUE), number)
  info <- getGenesFromList(input = names(list)[top_ids], filters = "entrezgene_id")
  info$diffusion_score <- list[as.character(info$entrezgene_id)]
  unique(info)
}

#' Degree-stratified test for a GO term
#'
#' @param go Character; GO term ID
#' @param network igraph object
#' @param attr Character; attribute in GO result to match (e.g. "entrezgene_id")
#' @param diffusion.scores Named numeric vector
#' @param perm Integer; number of permutations
#' @return Result of degreeStratifiedMonteCarlo()
#' @export
testGO <- function(go, network, attr = "entrezgene_id", diffusion.scores, perm = 1000) {
  go_df <- DownloadGenesInGO(request = go, attributes = attr)
  nodes <- igraph::V(network)$name %in% as.character(go_df[[attr]])
  degreeStratifiedMonteCarlo(network = network,
                             sample = igraph::V(network)$name[nodes],
                             diffusion.scores = diffusion.scores,
                             permutations = perm)
}

#' Get GO-annotated network nodes
#'
#' @inheritParams testGO
#' @return Character vector of vertex names in network
#' @export
getGONetworkNodes <- function(go, network, attr = "entrezgene_id") {
  go_df <- DownloadGenesInGO(request = go, attributes = attr)
  intersect(igraph::V(network)$name, as.character(go_df[[attr]]))
}

#' Intersection of top X percent lists
#'
#' @param network An igraph object.
#' @param list1 Named numeric vector
#' @param list2 Named numeric vector
#' @param percent Numeric; fraction of network to take
#' @return Character vector of intersecting external_gene_name
#' @import igraph
#' @export
topXpercentIntersection <- function(network, list1, list2, percent = 0.05) {
  cutoff <- as.integer(igraph::gorder(network) * percent)
  top1   <- getTopList(cutoff, list1)
  top2   <- getTopList(cutoff, list2)
  intersect(top1$external_gene_name, top2$external_gene_name)
}

#' Stratified sampling by a continuous control variable
#'
#' @param allNames Character vector of all IDs
#' @param sample Character vector of IDs to sample
#' @param controlVariable Numeric vector (same length as allNames)
#' @param controlSet Character vector; pool to sample from
#' @param sample_size Integer; size of sample
#' @param method Character; discretization method
#' @return Data.frame of sampled gene metadata
#' @import varrank
#' @export
stratifiedSample <- function(allNames, sample, controlVariable,
                             controlSet, sample_size, method = "sturges") {
  df <- logVariableDiscretisation(variable = controlVariable,
                                  names    = allNames,
                                  method   = method)
  ctrl_df <- subset(df, ID %in% controlSet)
  freq_tab <- merge(as.data.frame(table(subset(df, ID %in% sample)$LogVarDiscret)),
                    as.data.frame(table(df$LogVarDiscret)),
                    by = "Var1")
  freq_tab$Freq.fraction <- freq_tab$Freq.x / freq_tab$Freq.y
  df_prob <- merge(df, freq_tab[, c("Var1", "Freq.fraction")],
                   by.x = "LogVarDiscret", by.y = "Var1")
  sampled <- sample(df_prob$ID[df_prob$ID %in% controlSet],
                    size = sample_size,
                    prob = df_prob$Freq.fraction[df_prob$ID %in% controlSet])
  getGenesFromList(input = sampled, filters = "entrezgene_id")
}

#' Degree-stratified sampling for control selection
#'
#' @param network igraph object
#' @param sample Character vector of focal IDs
#' @param control Character vector of IDs to sample from
#' @param sample_size Integer; number to sample
#' @param method Character; discretization method
#' @return Data.frame of sampled gene metadata
#' @import igraph
#' @export
degreeStratifiedSample <- function(network, sample, control,
                                   sample_size, method = "sturges") {
  net_df  <- networkDiscretisation(network, method)
  ctrl_df <- subset(net_df, Gene.ID %in% control)
  freq_tab<- merge(as.data.frame(table(subset(net_df, Gene.ID %in% sample)$Category)),
                   as.data.frame(table(net_df$Category)),
                   by = "Var1")
  freq_tab$Freq.fraction <- freq_tab$Freq.x / freq_tab$Freq.y
  net_prob<- merge(net_df, freq_tab[, c("Var1", "Freq.fraction")],
                   by.x = "Category", by.y = "Var1")
  sampled <- sample(net_prob$Gene.ID[net_prob$Gene.ID %in% control],
                    size = sample_size,
                    prob = net_prob$Freq.fraction[net_prob$Gene.ID %in% control])
  getGenesFromList(input = sampled, filters = "entrezgene_id")
}
