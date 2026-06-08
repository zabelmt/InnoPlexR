# ============================================================
#  Structural analysis: clustering, components, Infomap export
# ============================================================

#' Local clustering coefficient for multiplex networks
#'
#' Computes the multiplex local clustering coefficient for each actor by
#' aggregating triangle counts across all layer pairs of the
#' supra-adjacency matrix using the block-tensor formulation.
#'
#' @param SupraAdjacencyMatrix A square sparse `dgCMatrix` as returned by
#'   [compute_supra_adjacency_matrix()].
#' @param Layers Integer. Number of layers.
#' @param Nodes Integer. Number of nodes per layer.
#'
#' @return A single-column numeric matrix of length `Nodes` with clustering
#'   coefficients in \eqn{[0, 1]}.
#'
#' @details
#' Requires `muxViz::SupraAdjacencyToBlockTensor()` to decompose the
#' \eqn{(NL \times NL)} supra-adjacency matrix into its
#' \eqn{L^2} blocks. Install muxViz with
#' `remotes::install_github("manlius/muxViz")` if not already available.
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' g1 <- igraph::sample_gnp(6, 0.6); igraph::V(g1)$name <- letters[1:6]
#' g2 <- igraph::sample_gnp(6, 0.6); igraph::V(g2)$name <- letters[1:6]
#' mx <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' S  <- compute_supra_adjacency_matrix(mx, omega = 1)
#' cl <- GetLocalClustering(S, Layers = 2, Nodes = 6)
#' }
GetLocalClustering <- function(SupraAdjacencyMatrix, Layers, Nodes) {

  if (!requireNamespace("muxViz", quietly = TRUE)) {
    stop("Package 'muxViz' is required for GetLocalClustering().\n",
         "Install it with: remotes::install_github(\"manlius/muxViz\")")
  }

  n <- Nodes * Layers

  # J - I: all-ones minus identity (potential triangle denominator)
  FMatrix <- Matrix::tcrossprod(
    Matrix::Matrix(1, n, 1, sparse = TRUE)
  ) - Matrix::Diagonal(n)

  # Numerator: closed triangles
  M3 <- SupraAdjacencyMatrix %*%
    SupraAdjacencyMatrix %*%
    SupraAdjacencyMatrix
  M3 <- methods::as(M3, "dgCMatrix")

  # Denominator: possible triangles excluding self-loops
  F3 <- SupraAdjacencyMatrix %*% FMatrix %*% SupraAdjacencyMatrix
  F3 <- methods::as(F3, "dgCMatrix")

  # Aggregate all L x L blocks into a single Nodes x Nodes matrix
  sum_blocks <- function(blocks) {
    Reduce("+", lapply(seq_len(Layers^2), function(k) blocks[[k]]))
  }

  B_num <- sum_blocks(muxViz::SupraAdjacencyToBlockTensor(M3, Layers, Nodes))
  B_den <- sum_blocks(muxViz::SupraAdjacencyToBlockTensor(F3, Layers, Nodes))

  idx  <- cbind(seq_len(Nodes), seq_len(Nodes))
  clus <- ifelse(B_den[idx] == 0, 0, B_num[idx] / B_den[idx])

  cbind(pmax(0, pmin(1, clus)))
}


#' Find connected components of a multiplex network
#'
#' Unions the intra-layer adjacency blocks of the supra-adjacency matrix,
#' binarises the result, and finds connected components via igraph.
#'
#' @param Multiplex A `Multiplex` object.
#' @param omega Non-negative numeric. Inter-layer coupling strength used
#'   when building the supra-adjacency matrix. Default `1`.
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{`no`}{Integer. Number of components.}
#'     \item{`csize`}{Integer vector of component sizes.}
#'     \item{`membership`}{Named integer vector mapping each node to its
#'       component ID.}
#'     \item{`graph`}{An igraph object on the `N`-node union graph.}
#'   }
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(8, 0.3); igraph::V(g1)$name <- paste0("n", 1:8)
#' g2 <- igraph::sample_gnp(8, 0.3); igraph::V(g2)$name <- paste0("n", 1:8)
#' mx <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' GetMultiComponents(mx)$no
GetMultiComponents <- function(Multiplex, omega = 1) {
  if (!is_Multiplex(Multiplex)) stop("Not a Multiplex object")

  N      <- Multiplex$Number_of_Nodes_Multiplex
  L      <- Multiplex$Number_of_Layers
  supra  <- compute_supra_adjacency_matrix(Multiplex, omega = omega)

  # Union of intra-layer blocks
  A_union <- Matrix::Matrix(0, nrow = N, ncol = N, sparse = TRUE)
  for (l in seq_len(L)) {
    idx     <- ((l - 1L) * N + 1L):(l * N)
    A_union <- A_union + supra[idx, idx]
  }
  A_union@x[] <- 1  # binarise

  g     <- igraph::graph_from_adjacency_matrix(A_union,
                                                mode     = "undirected",
                                                weighted = NULL,
                                                diag     = FALSE)
  comps <- igraph::components(g)

  list(
    no         = comps$no,
    csize      = comps$csize,
    membership = comps$membership,
    graph      = g
  )
}


#' Extract the largest connected component of a multiplex network
#'
#' A convenience wrapper around [GetMultiComponents()] that returns the
#' induced subgraph on the largest component's nodes.
#'
#' @param Multiplex A `Multiplex` object.
#' @param omega Non-negative numeric. Inter-layer coupling strength. Default
#'   `1`.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{`size`}{Integer. Number of nodes in the largest component.}
#'     \item{`nodes`}{Character vector of node names.}
#'     \item{`graph`}{An igraph object of the largest component.}
#'   }
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(8, 0.3); igraph::V(g1)$name <- paste0("n", 1:8)
#' g2 <- igraph::sample_gnp(8, 0.3); igraph::V(g2)$name <- paste0("n", 1:8)
#' mx  <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' lcc <- GetMultiLargestComponent(mx)
#' lcc$size
GetMultiLargestComponent <- function(Multiplex, omega = 1) {
  if (!is_Multiplex(Multiplex)) stop("Not a Multiplex object")

  comps      <- GetMultiComponents(Multiplex, omega = omega)
  largest_id <- which.max(comps$csize)
  nodes      <- names(which(comps$membership == largest_id))

  list(
    size  = comps$csize[largest_id],
    nodes = nodes,
    graph = igraph::induced_subgraph(comps$graph, nodes)
  )
}


#' Convert a Multiplex object to Infomap edge-list format
#'
#' Produces the inter- and intra-layer edge list required by the Infomap
#' community-detection algorithm, together with a node look-up table mapping
#' integer IDs back to original node names.
#'
#' @param x A `Multiplex` object.
#' @param omega Non-negative numeric. Inter-layer coupling weight. Default
#'   `1`.
#' @param weight_col Logical. Include a `weight` column in the output?
#'   Default `TRUE`.
#'
#' @return A list with two data frames:
#'   \describe{
#'     \item{`edge_list`}{Columns `layer_from`, `node_from`, `layer_to`,
#'       `node_to`, and (optionally) `weight`.}
#'     \item{`node_lookup`}{Columns `node_id` (integer) and `node_name`
#'       (character).}
#'   }
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(5, 0.5); igraph::V(g1)$name <- letters[1:5]
#' g2 <- igraph::sample_gnp(5, 0.5); igraph::V(g2)$name <- letters[1:5]
#' mx  <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' inf <- multiplex_to_infomap(mx)
#' head(inf$edge_list)
multiplex_to_infomap <- function(x, omega = 1, weight_col = TRUE) {

  if (!inherits(x, "Multiplex")) stop("Input must be a Multiplex object.")

  N            <- x$Number_of_Nodes_Multiplex
  L            <- x$Number_of_Layers
  P            <- x$Pool_of_Nodes
  Layers_Names <- names(x)[seq_len(L)]

  node_index <- stats::setNames(seq_len(N), P)

  # Intra-layer edges
  intra_list <- lapply(seq_len(L), function(k) {
    layer      <- x[[Layers_Names[k]]]
    has_weight <- "weight" %in% igraph::edge_attr_names(layer)
    el         <- igraph::as_edgelist(layer, names = TRUE)
    if (nrow(el) == 0L) return(NULL)

    df <- data.frame(
      layer_from = k,
      node_from  = node_index[el[, 1]],
      layer_to   = k,
      node_to    = node_index[el[, 2]],
      stringsAsFactors = FALSE
    )
    if (weight_col) df$weight <- if (has_weight) igraph::E(layer)$weight else 1
    df
  })

  # Inter-layer coupling edges
  inter_list <- NULL
  if (omega > 0 && L >= 2L) {
    layer_pairs <- utils::combn(seq_len(L), 2, simplify = FALSE)

    inter_list <- lapply(layer_pairs, function(pair) {
      k1 <- pair[1]; k2 <- pair[2]
      node_ids <- node_index[P]

      df <- data.frame(
        layer_from = c(rep(k1, N), rep(k2, N)),
        node_from  = rep(node_ids, 2),
        layer_to   = c(rep(k2, N), rep(k1, N)),
        node_to    = rep(node_ids, 2),
        stringsAsFactors = FALSE
      )
      if (weight_col) df$weight <- omega
      df
    })
  }

  all_links        <- do.call(rbind, c(intra_list, inter_list))
  rownames(all_links) <- NULL

  node_lookup <- data.frame(
    node_id   = seq_len(N),
    node_name = P,
    stringsAsFactors = FALSE
  )

  list(edge_list = all_links, node_lookup = node_lookup)
}
