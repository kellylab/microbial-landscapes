mapper.2.igraph <- function(mpr) {
  graf <- igraph::graph.adjacency(mpr$adjacency, mode="undirected")
  # name vertices so can track to subgraphs later
  igraph::V(graf)$name <- paste0("v", seq_along(igraph::V(graf)))
  graf
}