summary.TDAmapper <- function(m, plot = FALSE) {
  # points in each level
  pil <- data.frame(points = sapply(m$points_in_level, length))
  # points in each vertex
  ppv <- data.frame(points = sapply(m$points_in_vertex, length))
  # components
  graf <- igraph::graph.adjacency(m$adjacency, mode="undirected")
  csize <- igraph::components(graf)$csize
  if (plot) {
    p1 <- ggplot(pil, aes(x = points)) +
      geom_bar() +
      ggtitle("points in level")
    p2 <- ggplot(ppv, aes(x = points)) +
      geom_bar() +
      ggtitle("points in vertex")
    df <- data.frame(size = csize, rk = frank(-csize, ties.method = "first"))
    p3 <- ggplot(df, aes(x = rk, y = csize)) +
      geom_point() +
      labs(x = "rank", y = "size", title = "components")
    plot_grid(p1, p2, p3, ncol = 1, align = "v")

  } else {
    list(points.in.level = pil, points.per.vertex = ppv, component.size = csize)
  }
}