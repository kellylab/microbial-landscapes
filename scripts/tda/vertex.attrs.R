vertex.attrs <- function(graf, attr.df) {
  attr.names <- names(attr.df)
  for (name in attr.names) {
    graf <- igraph::set_vertex_attr(graf, name, V(graf), attr.df[[name]])
  }
  graf
}