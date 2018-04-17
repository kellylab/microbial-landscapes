#' Mapper vertex attribute
#'
#' Compute the attribute of a Mapper vertex by summarizing over the attributes
#' of its constituent data points.
#'
#' @param points a vector representing the data points within a vertex
#' @param point.attributes a vector of the point attributes
#' @param by.name index the point attribute vector by the names of the points
#' vector? If TRUE, both points and point.attributes must be named; if FALSE,
#' points must be numeric.
#' @param summ the function with which to summarize over point attributes
#'
#' @return the computed vertex attribute
#' @export
#'
#' @examples
vertex.attribute <- function(points, point.attributes,
                            by.name = TRUE, summ = "mean") {
  if (by.name) {
    v <- point.attributes[names(points)]
  } else {
    v <- point.attributes[points]
  }
  do.call(summ, list(x = v))
}