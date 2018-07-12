#' Plot the graph of a Mapper representation
#'
#' @param layout the Mapper graph as laid out with \code{ggraph}
#' @param node.aes Mapper vertex aesthetics, construct with e.g. \code{aes_}
#' @param labs list of plot
#' @param ... additional arguments to \code{geom_node_point}
#'
#' @return \code{ggraph} plot
#' @export
#'
#' @examples
plot.mapper <- function(layout, node.aes, labs = NULL, edge.color = "black",
                        ...) {
  p <- ggraph(layout)
  if (length(E(attr(layout, "graph"))) > 0) {
    p <- p + geom_edge_link0(colour = edge.color)
  }
  # ggraph(layout) +
    # geom_edge_link0() +
  p + geom_node_point(node.aes, ...) +
    labs(labs) +
    theme(aspect.ratio = 1) +
    theme_graph(base_family = "Helvetica")
}