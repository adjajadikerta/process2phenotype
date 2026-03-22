test_that("diffuseList runs on a toy graph with a named kernel", {
  library(igraph)

  g <- make_ring(5)
  V(g)$name <- as.character(1:5)

  K <- diag(5)
  rownames(K) <- V(g)$name
  colnames(K) <- V(g)$name

  scores <- diffuseList(
    network = g,
    inputGenes = "1",
    targetGenes = V(g)$name,
    kernel = K
  )

  expect_type(scores, "double")
  expect_equal(length(scores), 5)
  expect_equal(names(scores), V(g)$name)
  expect_equal(unname(scores), c(1, 0, 0, 0, 0))
})

test_that("diffuseList spreads signal to neighbours with a real kernel", {
  demo_net <- readRDS(
    system.file("extdata", "demo_network.rds", package = "process2phenotype")
  )
  demo_K <- readRDS(
    system.file("extdata", "demo_kernel.rds", package = "process2phenotype")
  )

  seed <- "5987"
  scores <- diffuseList(
    network    = demo_net,
    inputGenes = seed,
    kernel     = demo_K
  )

  expect_type(scores, "double")
  expect_equal(length(scores), igraph::vcount(demo_net))

  # Seed should have the highest score
  expect_equal(names(which.max(scores)), seed)

  # Direct neighbours should score higher on average than distant nodes
  nbrs      <- igraph::neighbors(demo_net, seed)$name
  non_nbrs  <- setdiff(names(scores), c(seed, nbrs))
  expect_gt(mean(scores[nbrs]), mean(scores[non_nbrs]))
})

test_that("diffuseList errors if input genes are not in target genes", {
  library(igraph)

  g <- make_ring(5)
  V(g)$name <- as.character(1:5)

  K <- diag(5)
  rownames(K) <- V(g)$name
  colnames(K) <- V(g)$name

  expect_error(
    diffuseList(
      network = g,
      inputGenes = "99",
      targetGenes = V(g)$name,
      kernel = K
    )
  )
})