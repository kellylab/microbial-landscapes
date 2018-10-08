local.extrema <- function(igraph_graf, min = TRUE, val = "knn") {
  is.extrema <- sapply(V(igraph_graf), function(v, graf) {
    x <- get.vertex.attribute(graf, val, v)
    neib.x <- get.vertex.attribute(graf, val, neighbors(graf, v))
    if (min) {
      all(x <= neib.x)
    } else {
      all(x >= neib.x)
    }
  }, graf = igraph_graf)
}