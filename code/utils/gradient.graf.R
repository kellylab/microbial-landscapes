gradient.graf <- function(in_graf, values, down = TRUE) {
  if (class(values) == "character") {
    values <- vertex_attr(in_graf, values)
  }
  ej <- as.data.frame(activate(in_graf, edges))
  from <- ej$from
  to <- ej$to
  delta <- values[to] - values[from]
  rectify <- function(from, to, delta, down) {
    newfrom <- from
    newto <- to
    if (down) {
      if (delta > 0) {
        newfrom <- to
        newto <- from
      }
    } else {
      if (delta < 0) {
        newfrom <- to
        newto <- from
      }
    }
    c(from = newfrom, to = newto, weight = abs(delta))
  }
  new.ej <- mapply(rectify, from = from, to = to, delta = delta,
                   MoreArgs = list(down = down))
  new.ej <- as.data.frame(t(new.ej))
  tbl_graph(nodes = as.data.frame(activate(in_graf, nodes)),
            edges = new.ej,
            directed = TRUE)
}