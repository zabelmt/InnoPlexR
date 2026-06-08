# InnoPlexR

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**InnoPlexR** provides tools for modelling innovation systems through multiplex
networks. It is designed for researchers integrating patent, publication, and
CORDIS funding data as separate network layers, enabling multilayer analysis of
knowledge flows, collaboration structures, and actor centrality across layers.
The package takes heavy inspiration from the muxViz package for multilayered
networks. Their work has been extended and adapted to work with large innovation
networks. 

If you use this package in academic work, please also consider citing the original
muxViz software and associated publications.

## Features

- `Multiplex` S3 class built on `igraph` layers
- Supra-adjacency, supra-Laplacian, and supra-transition matrix computation
- Centrality measures: eigenvector, PageRank, geodesic betweenness, random-walk betweenness
- Local clustering coefficient
- Connected component extraction
- Export to Infomap community-detection format

## Installation

Install from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("zabelmt/InnoPlexR")
```

### Dependencies

Core dependencies installed automatically:

- [`igraph`](https://igraph.org/r/) ≥ 1.3.0
- [`Matrix`](https://CRAN.R-project.org/package=Matrix) ≥ 1.5.0
- [`RSpectra`](https://CRAN.R-project.org/package=RSpectra) ≥ 0.16.0

Optional (for `GetLocalClustering()`):

```r
remotes::install_github("manlius/muxViz")
```

## Quick start

```r
library(InnoPlexR)

set.seed(42)
g_patents  <- igraph::sample_gnp(10, 0.35)
g_pubs     <- igraph::sample_gnp(10, 0.25)
g_cordis   <- igraph::sample_gnp(10, 0.30)

actor_names <- paste0("actor_", 1:10)
igraph::V(g_patents)$name <- igraph::V(g_pubs)$name <- igraph::V(g_cordis)$name <- actor_names

# Build multiplex object
mx <- create_Multiplex(g_patents, g_pubs, g_cordis,
                        Layers_Name = c("patents", "publications", "cordis"))
summary(mx)

# Supra-adjacency matrix
S <- compute_supra_adjacency_matrix(mx, omega = 1)

# Centrality measures
ec <- GetMultiEigenvectorCentrality(mx)
pr <- GetMultiPageRankCentrality(mx)
bc <- GetMultiGeodesicBetweenness(mx)

# Largest connected component
lcc <- GetMultiLargestComponent(mx)

# Export to Infomap format
inf <- multiplex_to_infomap(mx)
```

## Function reference

| Function | Description |
|---|---|
| `create_Multiplex()` | Build a `Multiplex` object from igraph layers |
| `is_Multiplex()` | Test if an object is a `Multiplex` |
| `compute_supra_adjacency_matrix()` | Supra-adjacency matrix |
| `compute_supra_laplacian_matrix()` | Supra-Laplacian matrix |
| `compute_supra_transition_matrix()` | Damped supra-transition matrix |
| `GetMultiEigenvectorCentrality()` | Multiplex eigenvector centrality |
| `GetMultiPageRankCentrality()` | Multiplex PageRank centrality |
| `GetMultiGeodesicBetweenness()` | Shortest-path betweenness |
| `GetMultiRandomWalkBetweenness()` | Random-walk betweenness (Monte Carlo) |
| `GetLocalClustering()` | Local clustering coefficient |
| `GetMultiComponents()` | Connected components |
| `GetMultiLargestComponent()` | Largest connected component |
| `multiplex_to_infomap()` | Export to Infomap edge-list format |

## References
De Domenico, M., Porter, M. A., & Arenas, A. (2015).
MuxViz: a tool for multilayer analysis and visualization of networks.

## License

MIT © Marcus Zabel
