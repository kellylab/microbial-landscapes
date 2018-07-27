gradient.graf <- function(in_graf, values, down = TRUE) {
  if (class(values) == "character") {
    values <- vertex_attr(in_graf, values)
  }
  ej.mat <- as.matrix(as.data.frame(activate(in_graf, edges)))
  ej.mat <- rbind(ej.mat, ej.mat[, c(2, 1)])
  keep <- apply(ej.mat, 1, function(vs, values, down) {
    if (down) {
      values[vs[1]] >= values[vs[2]]
    } else {
      values[vs[1]] <= values[vs[2]]
    }
  }, values = values, down = down)
  tbl_graph(nodes = as.data.frame(activate(in_graf, nodes)),
            edges = as.data.frame(ej.mat[keep, ]),
            directed = TRUE)
}