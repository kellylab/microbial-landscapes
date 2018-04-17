#' Create named vector
#'
#' Intended for use with vertex.attribute
#'
#' @param values vector of values
#' @param names vector of names
#'
#' @return a named vector
#' @export
#'
#' @examples
named.vector <- function(values, names) {
  v <- as.vector(values)
  names(v) <- names
  v
}