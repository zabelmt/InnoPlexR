# ============================================================
#  Multiplex S3 class — constructors, validators, helpers
# ============================================================

# ---- Internal helpers --------------------------------------------------

#' Add missing nodes to an igraph layer
#'
#' Ensures every node in the multiplex pool is present in a single layer
#' graph, inserting isolated vertices where needed.
#'
#' @param Layers An igraph object representing one layer.
#' @param Nr_Layers Integer. Total number of layers (unused internally but
#'   kept for interface compatibility).
#' @param NodeNames Character vector of all node names across all layers.
#'
#' @return The igraph object with any missing nodes added as isolates.
#' @keywords internal
add_missing_nodes <- function(Layers, Nr_Layers, NodeNames) {
  missing <- NodeNames[which(!NodeNames %in% igraph::V(Layers)$name)]
  if (length(missing) > 0) {
    Layers <- igraph::add_vertices(Layers, length(missing), name = missing)
  }
  Layers
}


# ---- Validators --------------------------------------------------------

#' Test whether an object is a Multiplex
#'
#' @param x Any R object.
#' @return `TRUE` if `x` has class `"Multiplex"`, `FALSE` otherwise.
#' @export
#' @examples
#' g <- igraph::make_ring(5)
#' igraph::V(g)$name <- letters[1:5]
#' mx <- create_Multiplex(g, Layers_Name = "ring")
#' is_Multiplex(mx)   # TRUE
#' is_Multiplex(g)    # FALSE
is_Multiplex <- function(x) {
  inherits(x, "Multiplex")
}


#' Test whether an object is a square sparse matrix
#'
#' Accepts both base `matrix` and `dgCMatrix` (sparse) objects.
#'
#' @param A A matrix object.
#' @return `TRUE` invisibly if the matrix is square; otherwise throws an
#'   error.
#' @export
#' @examples
#' A <- Matrix::Diagonal(3)
#' A <- as(A, "dgCMatrix")
#' is_square_dgCMatrix(A)   # TRUE
is_square_dgCMatrix <- function(A) {
  if (!inherits(A, "matrix") && !inherits(A, "dgCMatrix")) {
    stop("Supra_adjacency_Matrix must be a base matrix or sparse dgCMatrix.")
  }
  if (ncol(A) != nrow(A)) stop("Matrix must be square.")
  invisible(TRUE)
}


# ---- Constructor -------------------------------------------------------

#' Create a Multiplex network object
#'
#' Combines between one and six igraph layers into a single `Multiplex` S3
#' object. Layers that contain no edges are silently dropped. All node names
#' are pooled across every supplied layer to form the canonical node set.
#'
#' @param L1 An igraph object. **Required.** The first (or only) layer.
#' @param L2,L3,L4,L5,L6 Optional igraph objects for additional layers.
#'   Pass `NULL` (the default) to omit a layer.
#' @param Layers_Name Optional character vector of layer labels. Length must
#'   equal the number of non-empty layers. Defaults to
#'   `"Layer_1"`, `"Layer_2"`, …
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return A named list of class `"Multiplex"` containing:
#'   \describe{
#'     \item{`<layer name>`}{One igraph object per non-empty layer, with an
#'       edge attribute `type` set to the layer's name.}
#'     \item{`Pool_of_Nodes`}{Sorted character vector of all node names.}
#'     \item{`Number_of_Nodes_Multiplex`}{Integer. Total unique nodes.}
#'     \item{`Number_of_Layers`}{Integer. Number of non-empty layers.}
#'   }
#'
#' @export
#' @examples
#' set.seed(1)
#' g1 <- igraph::sample_gnp(10, 0.3)
#' g2 <- igraph::sample_gnp(10, 0.3)
#' igraph::V(g1)$name <- igraph::V(g2)$name <- paste0("n", 1:10)
#'
#' mx <- create_Multiplex(g1, g2, Layers_Name = c("A", "B"))
#' print(mx)
create_Multiplex <- function(L1, L2 = NULL, L3 = NULL, L4 = NULL,
                              L5 = NULL,  L6 = NULL,
                              Layers_Name, ...) {

  # --- Validate L1 -------------------------------------------------------
  if (!igraph::is_igraph(L1)) stop("Layer 1 is not an igraph object")

  # --- Replace NULLs with empty graphs -----------------------------------
  fill_empty <- function(L, label) {
    if (is.null(L)) return(igraph::make_empty_graph(n = 0, directed = FALSE))
    if (!igraph::is_igraph(L)) stop(label, " is not an igraph object")
    L
  }
  L2 <- fill_empty(L2, "Layer 2")
  L3 <- fill_empty(L3, "Layer 3")
  L4 <- fill_empty(L4, "Layer 4")
  L5 <- fill_empty(L5, "Layer 5")
  L6 <- fill_empty(L6, "Layer 6")

  # --- Keep only non-empty layers ----------------------------------------
  Layer_List <- Filter(function(g) igraph::ecount(g) >= 1,
                       list(L1, L2, L3, L4, L5, L6))

  # --- Pool of nodes (union of all vertex name attributes) ---------------
  all_names <- lapply(list(L1, L2, L3, L4, L5, L6),
                      function(g) igraph::V(g)$name)
  Pool_of_Nodes     <- sort(unique(unlist(all_names)))
  Number_of_Nodes   <- length(Pool_of_Nodes)
  Number_of_Layers  <- length(Layer_List)

  # --- Layer names -------------------------------------------------------
  if (missing(Layers_Name)) {
    Layers_Name <- paste0("Layer_", seq_len(Number_of_Layers))
  } else {
    if (!is.character(Layers_Name)) {
      stop("Layers_Name must be a character vector.")
    }
    if (length(Layers_Name) != Number_of_Layers) {
      stop("Length of Layers_Name (", length(Layers_Name),
           ") must equal the number of non-empty layers (", Number_of_Layers, ").")
    }
  }

  # --- Stamp edge attribute 'type' on each layer -------------------------
  Layer_List <- mapply(
    function(g, nm) igraph::set_edge_attr(g, "type", igraph::E(g), value = nm),
    Layer_List, Layers_Name,
    SIMPLIFY = FALSE
  )
  names(Layer_List) <- Layers_Name

  # --- Assemble and return -----------------------------------------------
  obj <- c(
    Layer_List,
    list(
      Pool_of_Nodes            = Pool_of_Nodes,
      Number_of_Nodes_Multiplex = Number_of_Nodes,
      Number_of_Layers         = Number_of_Layers
    )
  )
  class(obj) <- "Multiplex"
  obj
}


# ---- S3 methods --------------------------------------------------------

#' Print a Multiplex object
#'
#' @param x A `Multiplex` object.
#' @param ... Currently unused.
#' @return `x` invisibly.
#' @export
print.Multiplex <- function(x, ...) {
  cat("=== Multiplex network ===\n")
  cat("Layers :", x$Number_of_Layers, "\n")
  cat("Nodes  :", x$Number_of_Nodes_Multiplex, "\n\n")
  for (i in seq_len(x$Number_of_Layers)) {
    nm <- names(x)[i]
    cat("-- Layer", i, ":", nm, "--\n")
    print(x[[nm]])
    cat("\n")
  }
  invisible(x)
}

#' Summarise a Multiplex object
#'
#' @param object A `Multiplex` object.
#' @param ... Currently unused.
#' @return A data frame (invisibly) with one row per layer.
#' @export
summary.Multiplex <- function(object, ...) {
  layer_names <- names(object)[seq_len(object$Number_of_Layers)]
  df <- data.frame(
    layer   = layer_names,
    nodes   = vapply(layer_names, function(nm) igraph::vcount(object[[nm]]), integer(1)),
    edges   = vapply(layer_names, function(nm) igraph::ecount(object[[nm]]), integer(1)),
    stringsAsFactors = FALSE
  )
  cat("Multiplex network summary\n")
  cat("Total layers:", object$Number_of_Layers, "\n")
  cat("Total unique nodes:", object$Number_of_Nodes_Multiplex, "\n\n")
  print(df)
  invisible(df)
}
