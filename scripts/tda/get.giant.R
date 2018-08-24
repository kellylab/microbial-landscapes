get.giant <- function(tbl_graf) {
  library(tidygraph)
  library(igraph)
  comps <- components(tbl_graf)
  giant <- which(comps$csize == max(comps$csize))
  tbl_graf <- activate(tbl_graf, nodes)
  slice(tbl_graf, which(comps$membership == giant))
}
