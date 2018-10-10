vertex.2.points <- function(points.in.vertex) {
  out <- lapply(points.in.vertex, function(ps) {
    x <- data.table::data.table(point = ps)
    if (!is.null(names(ps))) {
      x$point.name <- names(ps)
    }
    x
    })
  out <- data.table::rbindlist(out, idcol = "vertex")
  out$vertex.name <- paste0("v", out$vertex)
  out
}