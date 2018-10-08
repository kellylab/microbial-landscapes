#' First k elements
#'
#' Sort a numeric vector and return the mean of the first k nonzero elements.
#' Useful for e.g. calculating k-nearest neighbor distance.
#'
#' @param v the vector
#' @param k how many elements to mean over
#' @param decreasing direction of sort
#'
#' @return the mean of the first k values
#' @export
#'
#' @examples
k.first <- function(v, k, decreasing = FALSE) {
  v <- sort(v[v > 0], decreasing = decreasing)
  mean(v[1:k])
}