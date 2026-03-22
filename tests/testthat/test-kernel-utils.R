# --- build_kernel tests -------------------------------------------------------

test_that("build_kernel returns a valid kernel from a small graph", {
  g <- igraph::make_ring(10)
  igraph::V(g)$name <- as.character(1:10)

  K <- build_kernel(g)

  expect_true(is.matrix(K))
  expect_true(is.numeric(K))
  expect_equal(nrow(K), 10)
  expect_equal(ncol(K), 10)
  expect_equal(rownames(K), igraph::V(g)$name)
  expect_equal(colnames(K), igraph::V(g)$name)
})

test_that("build_kernel passes ... to regularisedLaplacianKernel", {
  g <- igraph::make_ring(10)
  igraph::V(g)$name <- as.character(1:10)

  K_default    <- build_kernel(g)
  K_normalized <- build_kernel(g, normalized = TRUE)

  # Normalized and unnormalized kernels should differ

  expect_false(identical(K_default, K_normalized))
})

test_that("build_kernel errors on non-igraph input", {
  expect_error(build_kernel("not a graph"), "igraph")
})

test_that("build_kernel errors on unnamed vertices", {
  # make_ring without assigning names produces NULL V(g)$name
  g <- igraph::make_ring(5)
  expect_error(build_kernel(g), "names")
})

test_that("build_kernel errors on directed graph", {
  g <- igraph::make_ring(5, directed = TRUE)
  igraph::V(g)$name <- as.character(1:5)
  expect_error(build_kernel(g), "undirected")
})

test_that("build_kernel does not warn for small networks", {
  skip_on_cran()
  g <- igraph::make_ring(5)
  igraph::V(g)$name <- as.character(1:5)

  # Should not produce any warning about network size
  expect_no_warning(build_kernel(g), message = "nodes")
})

# --- load_kernel tests --------------------------------------------------------

test_that("load_kernel reads an .rds kernel", {
  demo_path <- system.file("extdata", "demo_kernel.rds",
                           package = "process2phenotype")
  K <- load_kernel(demo_path)

  expect_true(is.matrix(K))
  expect_true(is.numeric(K))
  expect_equal(nrow(K), ncol(K))
  expect_equal(rownames(K), colnames(K))
})

test_that("load_kernel warns about missing network nodes", {
  demo_path <- system.file("extdata", "demo_kernel.rds",
                           package = "process2phenotype")
  # Build a network that has an extra node not in the kernel
  demo_net <- readRDS(
    system.file("extdata", "demo_network.rds", package = "process2phenotype")
  )
  extra <- igraph::add_vertices(demo_net, 1, name = "FAKE_GENE")


  expect_warning(load_kernel(demo_path, network = extra), "not found in kernel")
})

test_that("load_kernel validates network argument is igraph", {
  demo_path <- system.file("extdata", "demo_kernel.rds",
                           package = "process2phenotype")
  expect_error(load_kernel(demo_path, network = "not a graph"), "igraph")
})

test_that("load_kernel reads an .rda kernel", {
  K_orig <- readRDS(
    system.file("extdata", "demo_kernel.rds", package = "process2phenotype")
  )
  tmp <- tempfile(fileext = ".rda")
  on.exit(unlink(tmp))
  save(K_orig, file = tmp)

  K <- load_kernel(tmp)
  expect_equal(K, K_orig)
})

test_that("load_kernel errors if .rda contains multiple objects", {
  a <- diag(3)
  rownames(a) <- colnames(a) <- letters[1:3]
  b <- a
  tmp <- tempfile(fileext = ".rda")
  on.exit(unlink(tmp))
  save(a, b, file = tmp)

  expect_error(load_kernel(tmp), "exactly 1 object")
})

test_that("load_kernel errors on non-matrix object", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(list(x = 1), tmp)

  expect_error(load_kernel(tmp), "not a numeric matrix")
})

test_that("load_kernel errors on non-square matrix", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  m <- matrix(1:6, nrow = 2, ncol = 3)
  rownames(m) <- c("a", "b")
  colnames(m) <- c("a", "b", "c")
  saveRDS(m, tmp)

  expect_error(load_kernel(tmp), "not square")
})

test_that("load_kernel errors on missing dimnames", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(diag(3), tmp)

  expect_error(load_kernel(tmp), "row names and column names")
})

test_that("load_kernel errors on mismatched dimnames", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  m <- diag(3)
  rownames(m) <- c("a", "b", "c")
  colnames(m) <- c("x", "y", "z")
  saveRDS(m, tmp)

  expect_error(load_kernel(tmp), "do not match")
})

test_that("load_kernel errors on missing file", {
  expect_error(load_kernel("/no/such/file.rds"), "not found")
})

test_that("load_kernel errors on unsupported extension", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines("x", tmp)

  expect_error(load_kernel(tmp), "Unsupported file extension")
})

# --- download_kernel tests ----------------------------------------------------

test_that("download_kernel errors on unknown name", {
  expect_error(download_kernel("nonexistent"), "Unknown kernel name")
})
