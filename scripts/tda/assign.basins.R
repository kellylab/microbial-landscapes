assign.basins <- function(tbl_graf, fn, down = TRUE) {
  library(tidygraph)
  library(igraph)
  source("local.extrema.R")
  # get local extrema according to specified function and direction
  tbl_graf <- mutate(tbl_graf, 
                     is.extremum = local.extrema(tbl_graf, val = fn, min = down))
  min2basin <- tbl_graf %>%  # assign adjacent minima to same basin
    activate(nodes) %>% 
    filter(is.extremum) %>%
    components %>%
    membership
  # convert to directed graph where edges go only in specified gradient direction
  source("gradient.graf.R") 
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
    mutate(basin = as.character(v2min[name]))
  tbl_graf
}