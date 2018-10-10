library(tidygraph)

merge.mpr.samples <- function(mpr, dt, by.y = "sample") {
  vertices <- mpr$graph %>%
    activate(nodes) %>%
    as.data.table %>%
    merge(mpr$mapping, by = "vertex") %>%
    merge(dt, by.x = "point", by.y = by.y)
  vertices
}