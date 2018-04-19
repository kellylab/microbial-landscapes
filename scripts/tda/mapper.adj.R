#' Fixed Mapper adjacency generation
#'
#' The adjacency matrix in \code{TDAmapper:mapper2d} only joins clusters in bin
#' \code{(i, j)} with those in bins \code{(i, j - 1)} and \code{(i - 1, j)}, but
#' not with those in bins \code{(i, j + 1)} and \code{(i + 1, j)}.
#' This function replaces that matrix.
#'
#' @param points.in.vertex \code{x$points_in_vertex}, where \code{x} is the
#' output of \code{TDAmapper::mapper2D}.
#'
#' @return a corrected adjacency matrix
#' @export
#'
#' @examples
mapper.adj <- function(points.in.vertex) {
  has.edge <- function(px, py) {
    length(intersect(px, py)) > 0
  }
  has.edge <- Vectorize(has.edge)
  outer(points.in.vertex, points.in.vertex, has.edge)
}
