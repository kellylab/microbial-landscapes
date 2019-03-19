#' Convert entries of vector to NA dependent on a Boolean vector
#'
#' @param bool a boolean vector
#' @param x a vector of same length as `bool`
#' @param invert if `TRUE`, invert values in `bool`
#'
#' @return
#' @export a copy of `x` with selected entries changed to `NA`
#'
#' @examples
nullify <- function(bool, x, invert = FALSE) {
  filler <- NA
  # class(filler) <- class(x)
  if (invert) bool <- !bool
  mapply(function(bi, xi, filler) if (bi) filler else xi,
         bi = bool, xi = x, MoreArgs = list(filler = filler))
}