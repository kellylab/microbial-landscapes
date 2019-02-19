library(igraph)
library(tidygraph)
library(ggraph)
#' Basic FR layout and plot of Mapper graph
#'
#' @param graf
#' @param edge
#' @param node
#' @param exclude.singletons
#' @param seed
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot.mapper.graph <- function(graf,
                              edge = geom_edge_link0(),
                              node = geom_node_point(),
                              exclude.singletons = FALSE,
                              seed = NULL,
                              layout = "fr",
                              ...) {
  set.seed(seed)
  if (exclude.singletons) {
    graf <- graf %>%
      activate(nodes) %>%
      filter(!in.singleton)
  }
  if (layout == "fr") {
    g <- ggraph(graf, layout, niter = 1000)
  } else {
    g <- ggraph(graf, layout, ...)
  }
  g +
    edge +
    node +
    theme_graph(base_family = "Helvetica") +
    theme(aspect.ratio = 1)
}