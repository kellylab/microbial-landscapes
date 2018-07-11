#' Layout `tbl_graph` with attached coordinates
#'
#' @param tbl_graph A `tbl_graph` with columns `x, y`.
#'
#' @return A `layout_igraph` `layout_ggraph` object
#' @export
#'
#' @examples
as.layout.manual <- function(tbl_graph) {
  tbl_graph <- tidygraph::activate(tbl_graph, nodes)
  d <- as.data.frame(tbl_graph)
  xy <- d[, c("x", "y")]
  ggraph::create_layout(tbl_graph, "manual", node.positions = xy)
}