#' @title Kernel utility functions
#' @description Helpers for building, loading, and downloading diffusion kernels.
#' @name kernel_utils
NULL

#' Build a diffusion kernel from an igraph network
#'
#' Thin wrapper around \code{\link[diffuStats]{regularisedLaplacianKernel}}.
#' Validates the input and warns if the network is large.
#'
#' @param network An undirected igraph object with named vertices.
#' @param ... Additional arguments passed to
#'   \code{\link[diffuStats]{regularisedLaplacianKernel}} (e.g. \code{sigma2},
#'   \code{add_diag}, \code{normalized}).
#' @return A square numeric kernel matrix with row and column names matching the
#'   vertex names of \code{network}.
#' @examples
#' library(igraph)
#' g <- make_ring(5)
#' V(g)$name <- as.character(1:5)
#' K <- build_kernel(g)
#' dim(K)
#' @export
build_kernel <- function(network, ...) {
  if (!inherits(network, "igraph")) {
    stop("'network' must be an igraph object")
  }
  if (is.null(igraph::V(network)$name)) {
    stop("'network' vertices must have names (V(network)$name)")
  }
  if (igraph::is_directed(network)) {
    stop("'network' must be undirected")
  }
  n <- igraph::vcount(network)
  if (n > 10000) {
    warning(
      "Network has ", n, " nodes. Kernel computation may require ",
      "substantial memory and time. See ?diffuStats::regularisedLaplacianKernel ",
      "for resource estimates.",
      call. = FALSE
    )
  }
  diffuStats::regularisedLaplacianKernel(graph = network, ...)
}

#' Load a precomputed diffusion kernel from file
#'
#' Reads a kernel matrix from an \code{.rds} or \code{.rda} file and validates
#' that it is a square numeric matrix with matching row and column names.
#'
#' @param path Character; path to an \code{.rds} or \code{.rda} file containing
#'   a kernel matrix.
#' @param network Optional igraph object. If supplied, a warning is issued for
#'   any network nodes absent from the kernel.
#' @return A square numeric kernel matrix.
#' @examples
#' \dontrun{
#' K <- load_kernel("path/to/kernel.rds")
#' K <- load_kernel("path/to/kernel.rds", network = my_graph)
#' }
#' @export
load_kernel <- function(path, network = NULL) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    K <- readRDS(path)
  } else if (ext %in% c("rda", "rdata")) {
    env <- new.env(parent = emptyenv())
    nms <- load(path, envir = env)
    if (length(nms) != 1L) {
      stop(
        "Expected .rda file to contain exactly 1 object, found ",
        length(nms), ": ", paste(nms, collapse = ", ")
      )
    }
    K <- env[[nms]]
  } else {
    stop("Unsupported file extension '.", ext, "'. Use .rds or .rda.")
  }

  if (!is.matrix(K) || !is.numeric(K)) {
    stop("Loaded object is not a numeric matrix")
  }
  if (nrow(K) != ncol(K)) {
    stop("Kernel matrix is not square (", nrow(K), " x ", ncol(K), ")")
  }
  if (is.null(rownames(K)) || is.null(colnames(K))) {
    stop("Kernel matrix must have both row names and column names")
  }
  if (!identical(rownames(K), colnames(K))) {
    stop("Kernel row names and column names do not match")
  }

  if (!is.null(network)) {
    if (!inherits(network, "igraph")) {
      stop("'network' must be an igraph object")
    }
    net_nodes <- igraph::V(network)$name
    missing <- setdiff(net_nodes, rownames(K))
    if (length(missing) > 0L) {
      warning(
        length(missing), " network node(s) not found in kernel: ",
        paste(utils::head(missing, 5), collapse = ", "),
        if (length(missing) > 5) paste0(", ... (", length(missing) - 5L, " more)"),
        call. = FALSE
      )
    }
  }
  K
}

#' Download a precomputed kernel from Zenodo
#'
#' Downloads a named kernel file from the project's Zenodo deposit and caches it
#' locally. If the file already exists and \code{overwrite} is \code{FALSE}, the
#' cached path is returned without re-downloading.
#'
#' @param name Character; one of \code{"human_biogrid"}, \code{"yeast_biogrid"},
#'   or \code{"biogrid_network"}.
#' @param dest_dir Character; directory in which to store the file. Defaults to
#'   \code{tools::R_user_dir("process2phenotype", which = "cache")}.
#' @param overwrite Logical; if \code{TRUE}, re-download even when a cached copy
#'   exists.
#' @return The file path (invisibly).
#' @examples
#' \dontrun{
#' path <- download_kernel("human_biogrid")
#' K <- load_kernel(path)
#' }
#' @export
download_kernel <- function(name, dest_dir = NULL, overwrite = FALSE) {
  file_map <- c(
    human_biogrid   = "BioGRID_Human_RLK.rds",
    yeast_biogrid   = "BioGRID_Yeast_RLK.rds",
    biogrid_network = "BioGRID_iGraph.rds"
  )
  if (!name %in% names(file_map)) {
    stop(
      "Unknown kernel name '", name, "'. ",
      "Choose one of: ", paste(names(file_map), collapse = ", ")
    )
  }
  filename <- file_map[[name]]
  base_url <- "https://zenodo.org/records/XXXXXXX/files"
  url <- paste0(base_url, "/", filename)

  if (is.null(dest_dir)) {
    dest_dir <- tools::R_user_dir("process2phenotype", which = "cache")
  }
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }
  dest <- file.path(dest_dir, filename)

  if (file.exists(dest) && !overwrite) {
    message("Using cached file: ", dest)
    return(invisible(dest))
  }

  message("Downloading ", filename, " from Zenodo...")
  utils::download.file(url, destfile = dest, mode = "wb")
  message("Saved to: ", dest)
  invisible(dest)
}
