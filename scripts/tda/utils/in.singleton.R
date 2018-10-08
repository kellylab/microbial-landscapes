#' Find Mapper components that represent only one data point
#'
#' @param point vector of oversampled data point indices or names
#' @param vertex vector of vertex mappings, same length as `point`
#' @param membership vector of component mappings, same length as
#' `unique(vertex)`
#'
#' @return Boolean vector of whether each vertex is in a singleton component
#' @export
#'
#' @examples
in.singleton <- function(point, vertex, membership) {
  df <- data.frame(point, vertex, mem = membership[vertex])
  df <- df %>%
    group_by(mem) %>%
    summarise(npts = length(unique(point))) %>%
    mutate(is.singleton = npts == 1)
  df$is.singleton[membership]
}