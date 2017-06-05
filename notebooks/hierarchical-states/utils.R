#' Generalized Jensen-Shannon divergence from matrix.
#'
#' @param x Matrix where each column is a discrete/binned probability distribution
#' (doesn't have to be normalized)
#' @param w weights of distributions, length must == ncol(x)
#'
#' @return JSD of all distribs in x
library(entropy)
# library(Matrix)
genJSD <- function(x, w = NULL) {
  # try to coerce x to a matrix
  if (!(class(x) %in% c("matrix", "dgeMatrix"))) {
    x <- as.matrix(x)
  }
  n <- ncol(x)
  if (is.null(w)) {
    w <- rep(1 / n, n)
  }
  # normalization
  x <- sweep(x, 2, colSums(x), FUN="/")
  # weighted sum of distributions
  wsum <- rowSums(x %*% w)
  # weighted entropy of distributions
  Hi <- sapply(1:n, function(i) {
    return(entropy.plugin(x[, i], unit = "log2"))
    })
  # the second term returns 1x1 matrix; as.numeric converts to scalar
  jsd <- entropy.plugin(wsum, unit = "log2") - as.numeric(Hi %*% w)
  return(jsd)
}

#' Given a circular ggraph layout, radially project the nodes outward.
#' Allows you to create radial layers or 'shells' of the same nodes formatted
#' different ways to show different data features.
#' Should you want to do such a crazy thing.
#'
#' @param layout a ggraph layout with circular = TRUE
#' @param dr the amount by which to project nodes outward
#'
#' @return the same layout with all nodes projected outward
#' @export
#'
#' @examples
RadialNudge <- function(layout, dr) {
  coord <- select(layout, x, y)
  coord <- mutate(coord, theta = atan2(y, x))
  coord <- mutate(coord, x = x + dr * cos(theta), y = y + dr * sin(theta))
  # dendrogram to igraph and back again
  layout <- create_layout(den_to_igraph(layout), "manual",
                          node.positions = coord
                          )
  return(layout)
}

#' Create ggraph layout using Rtsne
#'
#' @param graph the graph to be laid out (igraph)
#' @param dm distance matrix
#' @param seed seed the random number generator
#' @param perplexity
#'
#' @return layout_ggraph + layout_igraph
#' @export
#'
#' @examples
create_layout_tsne <- function(graph, dm, seed = NULL, perplexity = 15) {
  library(Rtsne)
  library(ggraph)
  library(igraph)
  set.seed(seed)
  ts <- Rtsne(dm, perplexity = perplexity)
  ts.layout <- create_layout(graph, "manual",
                             node.positions = data.table(x = ts$Y[, 1],
                                                         y = ts$Y[, 2]))
  return(ts.layout)
}

#' Makes an ggraph layout using non-metric (i.e. rank) multidimensional scaling
#' metaMDS from the vegan package.
#' Requires precomputed distance matrix.
#'
#' @param graph the graph to be laid out (igraph)
#' @param dm the distance matrix, must have vertex names in row and col names
#' @param ... additional arguments to metaMDS
#'
#' @return layout_ggraph + layout_igraph
#' @export
#'
#' @examples
create_layout_nmds <- function(graph, dm, ...) {
  library(vegan)
  library(data.table)
  nmds <- metaMDS(dm, ...)
  nmds <- data.table(x = nmds[["points"]][, "MDS1"], y = nmds[["points"]][, "MDS2"],
                     id = rownames(nmds[["points"]]))
  setkey(nmds, id)
  nmds <- nmds[V(graph)$name] # make sure coordinates in right order
  layout <- create_layout(graph, "manual", node.positions = nmds)
  return(layout)
}

#' Converts distance matrix to x, y coordinates using classical MDS.
#' Returns coordinates in a data.table and optionally eigenvalues.
#'
#' @param dm Distance matrix, must have row and column names.
#' @param eig Return eigenvalues?
#' @param ... Passed to `cmdscale`
#'
#' @return
#' If `eig`, list with point coordinates (`points`) and ranked eigenvalues (`eigrank`).
#' Otherwise just the coordinates. All as data.tables.
#'
#' @export
#'
#' @examples
CoordCMDS <- function(dm, eig = TRUE, ...) {
  fit <- cmdscale(dist(dm), eig = eig, ...)
  xform <- fit$points
  eigrank <- data.table(rank = seq(length(fit$eig)), value = fit$eig)
  colnames(xform) <- c("x", "y")
  rn <- rownames(xform)
  xform <- as.data.table(xform)
  xform[, sample := rn]
  if (eig) {
    list(points = xform, eigrank = eigrank)
  } else {
    xform
  }
}

#' Makes a named distance matrix from a data.table of (x, y, distance(x, y)).
#'
#' @param dt
#' @param id.x String, name of identifier column x
#' @param id.y String, name of identifier column y
#' @param value.cn name of column with distances
#'
#' @return
#' `dm`  Named distance `matrix`.
#'
#' @export
#'
#' @examples
MakeDistMatrix <- function(dt, id.x, id.y, value.cn = "distance") {
  formula <- paste(id.x, "~", id.y)
  dm <- dcast(dt, formula, value.var = value.cn)
  rn <- dm[, get(id.x)]
  dm <- as.matrix(dm[, -1])
  rownames(dm) <- rn
  dm
}

#' Given list of things, return data.table of non-redundant pairs between
#' things (including self-self).
#'
#' @param things vector of things
#' @param prefix string, prefix of column names (e.g. "thing")
#' @param suffix character length 2, suffix of column names (usually c("i", "j"), e.g.)
#' @param sep separator joining prefix and suffix (e.g. ".")
#'
#' @return data.table with thing i in one column and thing j in another
#' @export
#'
#' @examples
MakeNonRedundantPairs <- function(things, prefix = "thing",
                                  suffix = c("i", "j"), sep = ".") {
  library(data.table)
  N <- length(things)
  x <- lapply(1:N, function(n) {
    if (n == 1) {
      j <- things
    } else {
      j <- things[-seq(n - 1)]
    }
    out <- data.table(i = things[n], j = j)
    setnames(out, paste(prefix, suffix, sep = sep))
    out
  })
  x <- rbindlist(x)
  x
}