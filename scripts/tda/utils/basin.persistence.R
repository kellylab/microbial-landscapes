#' Persistence time/temporal correlation function of basins of attraction
#'
#' @param basin
#' @param time
#' @param scale If a time is associated with `n` basins, set the contribution of
#' each to `1/n`
#'
#' @return
#' @export
#'
#' @examples
basin.persistence <- function(basin, time, scale = FALSE, summarize = FALSE) {
  library(tidyverse)
  df <- data.frame(basin = basin, time = time, foo = 1)
  df <- unique(df)
  m <- reshape2::acast(df, time ~ basin, value.var = "foo", fill = 0)
  if (scale) {
    m <- sweep(m, 1, rowSums(m), "/")
  }
  time <- as.numeric(rownames(m))
  cor.fun <- function(has.basin, time, summarize = FALSE) {
    names(time) <- time
    time.pairs <- outer(time, time, "-")
    time.pairs <- reshape2::melt(time.pairs, varnames = c("time.j", "time.i"),
                       value.name = "delta.t")
    time.pairs <- dplyr::filter(time.pairs, delta.t > 0)
    names(has.basin) <- time
    time.pairs <- dplyr::mutate(time.pairs,
                                has.i = has.basin[as.character(time.i)],
                                has.j = has.basin[as.character(time.j)])
    time.pairs <- dplyr::filter(time.pairs, has.i > 0)
    if (summarize) {
      time.pairs <- dplyr::group_by(time.pairs, delta.t)
      cf <- dplyr::summarise(time.pairs, f = mean(has.j), var = var(has.j),
                             n = n())
    } else {
      cf <- dplyr::select(time.pairs, delta.t, f = has.j)
    }
    cf
  }
  result <- apply(m, 2, cor.fun, time = time, summarize = summarize)
  dplyr::bind_rows(result, .id = "basin")
}