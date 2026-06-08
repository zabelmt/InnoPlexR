#' InnoPlexR: Build and Analyze Multiplex Networks
#'
#' @description
#' Tools for constructing and analyzing multiplex networks. The package
#' provides a `Multiplex` S3 class built on top of igraph layers, along with
#' functions for computing matrix representations (supra-adjacency, supra-
#' Laplacian, supra-transition) and a suite of multilayer centrality measures.
#'
#' @section Core workflow:
#' 1. Build a `Multiplex` object with [create_Multiplex()].
#' 2. Compute matrix representations with [compute_supra_adjacency_matrix()],
#'    [compute_supra_laplacian_matrix()], or [compute_supra_transition_matrix()].
#' 3. Compute centralities: [GetMultiEigenvectorCentrality()],
#'    [GetMultiPageRankCentrality()], [GetMultiGeodesicBetweenness()],
#'    [GetMultiRandomWalkBetweenness()].
#' 4. Analyse structure: [GetMultiComponents()], [GetMultiLargestComponent()],
#'    [GetLocalClustering()].
#' 5. Export for external tools: [multiplex_to_infomap()].
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom igraph is_igraph make_empty_graph add_vertices V E ecount
#'   set_edge_attr edge_attr_names as_adjacency_matrix as_edgelist
#'   graph_from_adjacency_matrix betweenness components induced_subgraph
#' @importFrom Matrix Diagonal bdiag kronecker Matrix rowSums colSums tcrossprod
#' @importFrom RSpectra eigs_sym
#' @importFrom methods as is new
#' @importFrom utils combn head tail
## usethis namespace: end
NULL
