# ============================================================
#  Supra-matrix representations
# ============================================================

#' Compute the supra-adjacency matrix of a multiplex network
#'
#' Builds the block \eqn{(N \cdot L) \times (N \cdot L)} supra-adjacency
#' matrix by stacking the per-layer adjacency matrices on the block diagonal
#' and filling the off-diagonal coupling blocks with the scalar `omega`.
#'
#' @param x A `Multiplex` object created by [create_Multiplex()].
#' @param omega Non-negative numeric. Inter-layer coupling strength. Default
#'   `167` (a strong coupling that keeps components connected).
#'
#' @return A sparse `dgCMatrix` of dimension \eqn{(N \cdot L) \times
#'   (N \cdot L)}, with row/column names of the form `"nodeName_layerIndex"`.
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(5, 0.5); igraph::V(g1)$name <- letters[1:5]
#' g2 <- igraph::sample_gnp(5, 0.5); igraph::V(g2)$name <- letters[1:5]
#' mx <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' S  <- compute_supra_adjacency_matrix(mx, omega = 1)
#' dim(S)
compute_supra_adjacency_matrix <- function(x, omega = 167) {

  if (!is_Multiplex(x)) stop("Not a Multiplex object")
  if (!is.numeric(omega) || length(omega) != 1 || omega < 0) {
    stop("omega must be a non-negative numeric scalar.")
  }

  N            <- x$Number_of_Nodes_Multiplex
  L            <- x$Number_of_Layers
  P            <- x$Pool_of_Nodes
  Layers_Names <- names(x)[seq_len(L)]

  if (L < 2) stop("At least 2 layers are required.")

  # Build per-layer N x N adjacency matrices aligned to Pool_of_Nodes order
  Layers_List <- lapply(Layers_Names, function(layer_name) {
    layer <- x[[layer_name]]

    if ("weight" %in% igraph::edge_attr_names(layer)) {
      Adj <- igraph::as_adjacency_matrix(layer, sparse = TRUE, attr = "weight")
    } else {
      Adj <- igraph::as_adjacency_matrix(layer, sparse = TRUE)
    }

    layer_nodes   <- rownames(Adj)
    missing_nodes <- setdiff(P, layer_nodes)

    if (length(missing_nodes) > 0) {
      n_miss   <- length(missing_nodes)
      zero_col <- Matrix::Matrix(0, nrow = nrow(Adj), ncol = n_miss, sparse = TRUE)
      zero_row <- Matrix::Matrix(0, nrow = n_miss,    ncol = N,      sparse = TRUE)
      colnames(zero_col) <- missing_nodes
      rownames(zero_row) <- missing_nodes
      Adj <- cbind(Adj, zero_col)
      Adj <- rbind(Adj, zero_row)
    }

    methods::as(Adj[P, P], "dgCMatrix")
  })

  # Block-diagonal intra-layer part
  B <- Matrix::bdiag(Layers_List)

  # Off-diagonal inter-layer coupling
  Cmat        <- matrix(omega, nrow = L, ncol = L)
  diag(Cmat)  <- 0
  Inter       <- Matrix::kronecker(Cmat, Matrix::Diagonal(N))

  Supra <- B + Inter

  # Row/col labels: nodeName_layerIndex
  vnames      <- rownames(Layers_List[[1]])
  if (is.null(vnames)) vnames <- as.character(seq_len(N))
  block_names <- unlist(lapply(seq_len(L), function(k) paste0(vnames, "_", k)))
  dimnames(Supra) <- list(block_names, block_names)

  methods::as(Supra, "dgCMatrix")
}


#' Compute the supra-Laplacian matrix of a multiplex network
#'
#' Returns \eqn{L = D - A} where \eqn{A} is the supra-adjacency matrix
#' (see [compute_supra_adjacency_matrix()]) and \eqn{D} is the diagonal
#' degree matrix of the supra-graph.
#'
#' @param Multiplex A `Multiplex` object.
#' @param omega Non-negative numeric. Inter-layer coupling strength passed to
#'   [compute_supra_adjacency_matrix()]. Default `167`.
#'
#' @return A sparse `dgCMatrix` Laplacian of dimension
#'   \eqn{(N \cdot L) \times (N \cdot L)}.
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(5, 0.5); igraph::V(g1)$name <- letters[1:5]
#' g2 <- igraph::sample_gnp(5, 0.5); igraph::V(g2)$name <- letters[1:5]
#' mx <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' L  <- compute_supra_laplacian_matrix(mx)
compute_supra_laplacian_matrix <- function(Multiplex, omega = 167) {
  supra_A <- compute_supra_adjacency_matrix(Multiplex, omega = omega)
  degree  <- Matrix::rowSums(supra_A)
  D       <- Matrix::Diagonal(x = degree)
  D - supra_A
}


#' Compute the supra-transition (PageRank) matrix of a multiplex network
#'
#' Column-normalises the supra-adjacency matrix and applies the standard
#' PageRank damping formula:
#' \deqn{T = d \cdot \hat{A} + (1-d)/N \cdot \mathbf{1}\mathbf{1}^T}
#' where \eqn{\hat{A}} is the column-stochastic transition matrix and
#' \eqn{d} is the damping factor.
#'
#' @param Multiplex A `Multiplex` object.
#' @param damping Numeric in \eqn{(0,1)}. PageRank damping factor. Default
#'   `0.85`.
#' @param omega Non-negative numeric. Inter-layer coupling strength. Default
#'   `1`.
#'
#' @return A dense numeric matrix of dimension
#'   \eqn{(N \cdot L) \times (N \cdot L)}.
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(5, 0.5); igraph::V(g1)$name <- letters[1:5]
#' g2 <- igraph::sample_gnp(5, 0.5); igraph::V(g2)$name <- letters[1:5]
#' mx <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' T_mat <- compute_supra_transition_matrix(mx)
compute_supra_transition_matrix <- function(Multiplex, damping = 0.85, omega = 1) {

  n_nodes  <- Multiplex$Number_of_Nodes_Multiplex
  n_layers <- Multiplex$Number_of_Layers
  N        <- n_nodes * n_layers

  supra_A   <- compute_supra_adjacency_matrix(Multiplex, omega = omega)
  col_sums  <- Matrix::colSums(supra_A)
  col_sums[col_sums == 0] <- 1          # handle dangling nodes

  T_matrix <- Matrix::t(Matrix::t(supra_A) / col_sums)
  teleport_matrix <- Matrix::Matrix(1 / N, nrow = N, ncol = N, sparse = FALSE)

  damping * T_matrix + (1 - damping) * teleport_matrix
}
