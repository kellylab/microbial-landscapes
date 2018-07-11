#' Attach layout coordinates to `tbl_graph`
#'
#' @param tbl_graph A `tbl_graph` object.
#' @param ... Arguments to `ggraph::create_layout`
#'
#' @return A `tbl_graph` with columns `x, y`
#' @export
#'
#' @examples
layout.tbl.graph <- function(tbl_graph, ...) {
  lo <- ggraph::create_layout(tbl_graph, ...)
  df <- as.data.frame(lo)
  tbl_graph <- tidygraph::activate(tbl_graph, nodes)
  tidygraph::mutate(tbl_graph, x = lo$x, y = lo$y)
}