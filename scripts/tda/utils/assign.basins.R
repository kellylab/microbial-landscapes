#' Find local minima and basins of attractions in Mapper graph
#'
#' @param tbl_graf Mapper graph as tidygraph::tbl_graph object
#' @param fn vertex attribute representing value of quasipotential function
#' @param down minima and descent if TRUE; maxima and ascent if FALSE
#' @param ignore.singletons set basins with onl 1 vertex to NA
#' @param giant.only only find basins in giant component, set others to NA
#'
#' @return updated tbl_graph object
#' @export
#'
#' @examples
assign.basins <- function(tbl_graf, fn, down = TRUE, ignore.singletons = TRUE,
                          giant.only = FALSE) {
  library(tidygraph)
  library(igraph)
  # source("local.extrema.R")
  # source("gradient.graf.R")
  # source("get.giant.R")
  orig.graf <- tbl_graf
  if (giant.only) {
    tbl_graf <- get.giant(tbl_graf)
  }
  # get local extrema according to specified function and direction
  tbl_graf <- mutate(tbl_graf,
                     is.extremum = local.extrema(tbl_graf, val = fn, min = down))
  min2basin <- tbl_graf %>%  # assign adjacent minima to same basin
    activate(nodes) %>%
    filter(is.extremum) %>%
    components %>%
    membership
  # convert to directed graph where edges go only in specified gradient direction
  dgraf <- gradient.graf(tbl_graf, fn, down = down)
  # directed graph distance of each vertex to each extremum
  dist.2.min <- distances(dgraf, to = which(V(tbl_graf)$is.extremum), weights = NA)
  # function to map each vertex to a basin from distances
  get.basin <- function(v, lab = names(v), map) {
    w <- lab[v == min(v)]
    b <- unique(map[w])
    if (length(b) == 1) {
      b
    } else {
      NA
    }
  }
  v2min <- apply(dist.2.min, 1, get.basin, lab = colnames(dist.2.min),
                 map = min2basin)
  tbl_graf <- tbl_graf %>%
    activate(nodes) %>%
    mutate(basin = v2min[name])
  if (ignore.singletons) {
    tbl_graf <- tbl_graf %>%
      group_by(basin) %>%
      mutate(is.extremum = if (n() > 1) is.extremum else FALSE,
             newbasin = if (n() > 1) basin else NA) %>%
      ungroup %>%
      mutate(basin = newbasin) %>%
      mutate(newbasin = NULL)
  }
  if (giant.only) {
    df <- as.data.frame(select(tbl_graf, vertex, is.extremum, basin))
    tbl_graf <- left_join(orig.graf, df, by = "vertex")
  }
  tbl_graf %>%
    activate(nodes) %>%
    mutate(basin, as.factor(basin))
  tbl_graf
}