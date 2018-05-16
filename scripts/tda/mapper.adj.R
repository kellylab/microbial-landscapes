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
  nv <- length(points.in.vertex)
  np <- max(sapply(points.in.vertex, max)) # max point id
  vp <- matrix(0, nv, np)
  for (vi in seq_len(nv)) {
    vp[vi, points.in.vertex[[vi]]] <- 1
  }
  (vp %*% t(vp)) > 0
}
