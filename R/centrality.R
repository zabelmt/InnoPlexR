# ============================================================
#  Multiplex centrality measures
# ============================================================

#' Multiplex eigenvector centrality
#'
#' Computes the leading eigenvector of the supra-adjacency matrix and maps
#' the resulting scores back to actors across layers.
#'
#' @param Multiplex A `Multiplex` object created by [create_Multiplex()].
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{`by_layer`}{Numeric matrix (`nodes x layers`) of per-layer
#'       eigenvector centrality scores.}
#'     \item{`aggregate`}{Named numeric vector of layer-summed scores
#'       normalised to \eqn{[0, 1]}.}
#'     \item{`eigenvalue`}{The leading eigenvalue.}
#'   }
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(8, 0.4); igraph::V(g1)$name <- paste0("n", 1:8)
#' g2 <- igraph::sample_gnp(8, 0.4); igraph::V(g2)$name <- paste0("n", 1:8)
#' mx <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' ec <- GetMultiEigenvectorCentrality(mx)
#' head(ec$aggregate)
GetMultiEigenvectorCentrality <- function(Multiplex) {

  n_nodes  <- Multiplex$Number_of_Nodes_Multiplex
  n_layers <- Multiplex$Number_of_Layers

  supra_A  <- compute_supra_adjacency_matrix(Multiplex)
  eig      <- RSpectra::eigs_sym(supra_A, k = 1, which = "LM")
  ec       <- abs(as.vector(eig$vectors))

  ec_matrix <- matrix(ec,
                      nrow = n_nodes, ncol = n_layers,
                      dimnames = list(
                        rownames(supra_A)[seq_len(n_nodes)],
                        paste0("layer_", seq_len(n_layers))
                      ))

  ec_agg  <- rowSums(ec_matrix)
  list(
    by_layer   = ec_matrix,
    aggregate  = ec_agg / max(ec_agg),
    eigenvalue = eig$values
  )
}


#' Multiplex PageRank centrality
#'
#' Runs power iteration on the damped supra-transition matrix to compute
#' multilayer PageRank scores.
#'
#' @param Multiplex A `Multiplex` object.
#' @param damping Numeric in \eqn{(0,1)}. PageRank damping factor. Default
#'   `0.85`.
#' @param omega Non-negative numeric. Inter-layer coupling strength. Default
#'   `1`.
#' @param tol Convergence tolerance for power iteration. Default `1e-6`.
#' @param max_iter Maximum number of power-iteration steps. Default `100`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`by_layer`}{Numeric matrix (`nodes x layers`).}
#'     \item{`aggregate`}{Named numeric vector normalised to \eqn{[0, 1]}.}
#'   }
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(8, 0.4); igraph::V(g1)$name <- paste0("n", 1:8)
#' g2 <- igraph::sample_gnp(8, 0.4); igraph::V(g2)$name <- paste0("n", 1:8)
#' mx <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' pr <- GetMultiPageRankCentrality(mx)
#' head(pr$aggregate)
GetMultiPageRankCentrality <- function(Multiplex, damping = 0.85, omega = 1,
                                        tol = 1e-6, max_iter = 100) {

  n_nodes  <- Multiplex$Number_of_Nodes_Multiplex
  n_layers <- Multiplex$Number_of_Layers
  N        <- n_nodes * n_layers

  T_damped <- compute_supra_transition_matrix(Multiplex,
                                               damping = damping,
                                               omega   = omega)

  pr <- rep(1 / N, N)

  for (iter in seq_len(max_iter)) {
    pr_new <- as.vector(T_damped %*% pr)
    pr_new <- pr_new / sum(pr_new)

    if (max(abs(pr_new - pr)) < tol) {
      message("PageRank converged in ", iter, " iterations.")
      break
    }
    pr <- pr_new
  }

  pr_matrix <- matrix(pr,
                      nrow = n_nodes, ncol = n_layers,
                      dimnames = list(
                        rownames(T_damped)[seq_len(n_nodes)],
                        paste0("layer_", seq_len(n_layers))
                      ))

  pr_agg <- rowSums(pr_matrix)
  list(
    by_layer  = pr_matrix,
    aggregate = pr_agg / max(pr_agg)
  )
}


#' Multiplex geodesic betweenness centrality
#'
#' Computes shortest-path betweenness on the supra-graph using
#' [igraph::betweenness()]. Edge weights are inverted so that stronger
#' connections correspond to shorter distances.
#'
#' @param Multiplex A `Multiplex` object.
#' @param normalized Logical. Normalise betweenness scores? Default `TRUE`.
#' @param weighted Logical. Use edge weights (inverted)? Default `TRUE`.
#' @param omega Non-negative numeric. Inter-layer coupling strength. Default
#'   `1`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`by_layer`}{Numeric matrix (`nodes x layers`).}
#'     \item{`aggregate`}{Named numeric vector normalised to \eqn{[0, 1]}.}
#'   }
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(6, 0.5); igraph::V(g1)$name <- letters[1:6]
#' g2 <- igraph::sample_gnp(6, 0.5); igraph::V(g2)$name <- letters[1:6]
#' mx <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' bc <- GetMultiGeodesicBetweenness(mx)
#' head(bc$aggregate)
GetMultiGeodesicBetweenness <- function(Multiplex, normalized = TRUE,
                                         weighted = TRUE, omega = 1) {

  n_nodes  <- Multiplex$Number_of_Nodes_Multiplex
  n_layers <- Multiplex$Number_of_Layers

  supra_A <- compute_supra_adjacency_matrix(Multiplex, omega)

  g <- igraph::graph_from_adjacency_matrix(
    supra_A,
    mode     = "undirected",
    weighted = if (weighted) TRUE else NULL,
    diag     = FALSE
  )

  if (weighted) igraph::E(g)$weight <- 1 / igraph::E(g)$weight

  bc_raw <- igraph::betweenness(
    g,
    v          = igraph::V(g),
    directed   = FALSE,
    weights    = if (weighted) igraph::E(g)$weight else NA,
    normalized = normalized
  )

  bc_matrix <- matrix(bc_raw,
                      nrow = n_nodes, ncol = n_layers,
                      dimnames = list(
                        colnames(supra_A)[seq_len(n_nodes)],
                        paste0("layer_", seq_len(n_layers))
                      ))

  bc_agg <- rowSums(bc_matrix)
  list(
    by_layer  = bc_matrix,
    aggregate = bc_agg / max(bc_agg)
  )
}


#' Multiplex random-walk betweenness centrality (Monte Carlo)
#'
#' Approximates current-flow (random-walk) betweenness via a Monte Carlo
#' approach: the pseudoinverse of the supra-Laplacian is computed once, then
#' node voltages are solved for a random sample of source–target pairs. This
#' is far more scalable than the exact \eqn{O(N^3)} computation.
#'
#' @param Multiplex A `Multiplex` object.
#' @param tol Numeric. Eigenvalue threshold for truncating the pseudoinverse.
#'   Default `1e-6`.
#' @param normalized Logical. Reserved for interface consistency; currently
#'   scores are always normalised to \eqn{[0, 1]}.
#' @param n_samples Integer. Number of random source–target pairs to sample.
#'   Larger values give more accurate estimates at higher cost. Default
#'   `1000`.
#' @param seed Integer. Random seed for reproducibility. Default `42`.
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{`by_layer`}{Numeric matrix (`nodes x layers`) of raw scores.}
#'     \item{`aggregate`}{Named numeric vector normalised to \eqn{[0, 1]}.}
#'     \item{`n_samples`}{The number of pairs that were sampled.}
#'     \item{`total_pairs`}{Total possible pairs \eqn{N(N-1)/2} (for context).}
#'   }
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' g1 <- igraph::sample_gnp(8, 0.4); igraph::V(g1)$name <- paste0("n", 1:8)
#' g2 <- igraph::sample_gnp(8, 0.4); igraph::V(g2)$name <- paste0("n", 1:8)
#' mx  <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' rwb <- GetMultiRandomWalkBetweenness(mx, n_samples = 100)
#' head(rwb$aggregate)
#' }
GetMultiRandomWalkBetweenness <- function(Multiplex, tol = 1e-6,
                                           normalized = TRUE,
                                           n_samples = 1000,
                                           seed = 42) {

  n_nodes  <- Multiplex$Number_of_Nodes_Multiplex
  n_layers <- Multiplex$Number_of_Layers
  N        <- n_nodes * n_layers

  supra_A <- compute_supra_adjacency_matrix(Multiplex)
  L_mat   <- compute_supra_laplacian_matrix(Multiplex)

  # Pseudoinverse via truncated eigendecomposition
  k        <- min(N - 1L, 100L)
  eig      <- RSpectra::eigs_sym(L_mat, k = k, which = "LM")
  non_zero <- abs(eig$values) > tol
  eig_vals <- eig$values[non_zero]
  eig_vecs <- eig$vectors[, non_zero, drop = FALSE]
  L_pinv   <- eig_vecs %*% Matrix::Diagonal(x = 1 / eig_vals) %*% t(eig_vecs)

  # Sample random source-target pairs
  set.seed(seed)
  pairs <- replicate(n_samples, sample.int(N, 2L), simplify = FALSE)

  rwb <- numeric(N)

  for (pair in pairs) {
    s <- pair[1]; t <- pair[2]

    b_st     <- numeric(N)
    b_st[s]  <-  1
    b_st[t]  <- -1

    voltages <- as.vector(L_pinv %*% b_st)

    for (v in seq_len(N)) {
      nbrs       <- which(supra_A[v, ] > 0)
      edge_flows <- supra_A[v, nbrs] * (voltages[v] - voltages[nbrs])
      rwb[v]     <- rwb[v] + sum(abs(edge_flows)) / 2
    }
  }

  # Scale up to full pair count
  total_pairs <- N * (N - 1L) / 2L
  rwb         <- rwb * (total_pairs / n_samples)

  rwb_matrix <- matrix(rwb,
                       nrow = n_nodes, ncol = n_layers,
                       dimnames = list(
                         rownames(supra_A)[seq_len(n_nodes)],
                         paste0("layer_", seq_len(n_layers))
                       ))

  rwb_agg <- rowSums(rwb_matrix)
  list(
    by_layer    = rwb_matrix,
    aggregate   = rwb_agg / max(rwb_agg),
    n_samples   = n_samples,
    total_pairs = total_pairs
  )
}
