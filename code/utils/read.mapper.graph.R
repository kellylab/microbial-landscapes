library(data.table)
library(tidygraph)

read.mapper.graph <- function(directory) {
  vertices <- fread(paste0(directory, "/vertices.txt"))
  edges <- fread(paste0(directory, "/edges.txt"))
  graf <- tbl_graph(vertices, edges, FALSE)
  v2p <- fread(paste0(directory, "/vertices-to-points.txt"))
  list(graph = graf, mapping = v2p)
}