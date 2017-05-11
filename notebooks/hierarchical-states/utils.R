NormDegrees <- function(matrix) {
  # prevent bugging out on 1x1 mats
  matrix <- as.matrix(matrix)
  nhosts <- dim(matrix)[1]
  nphages <- dim(matrix)[2]
  host.degree.norm <- rowSums(matrix) / nphages
  phage.degree.norm <- colSums(matrix) / nhosts
  return(list('host'=host.degree.norm, 'phage'=phage.degree.norm))
}

Connectivity <- function(matrix) {
  matrix <- as.matrix(matrix)
  md <- dim(matrix)
  if (md[1] == 0 | md[2] ==0) {
    c <- 0
  }
  else {
    c <- sum(matrix) / (md[1] * md[2])
  }
  return(c)
}

# get the size fraction (um) from which the host was collected by ID
HostSizeFraction <- function(x) {
  id <- unlist(strsplit(x[["member"]], "[.]"))[3]
  if (id %in% c(54, 55, 56)) {
    sf <- "63"
  }
  else if (id %in% c(51, 52, 53)) {
    sf <- "5"
  }
  else if (id %in% c(48, 49, 50)) {
    sf <- "1"
  }
  else if (id %in% c(45, 46, 47)) {
    sf <- "0.2"
  }
  else {
    sf <- NA
  }
  return(sf)
}

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

#' Function to output phage and protein eigenvector centralities.
#' Protein centrality calculated from phage centrality.
#' See Faust (1997), MEJ Newman (2004).
#'
#' @param u.scores dataframe/table with columns phage, host.prot, and u.phage
#' @param norm normalize centrality scores? (currently unimplemented)
#'
#' @return list of 2 data frames: (phage centrality, protein centrality)
#' @export
#'
#' @examples
PPCentrality <- function(u.scores, by = "host.prot", norm=FALSE) {
  # clear unnecessary columns and make a new data.table
  u.work <- select(u.scores, phage, get(by), u.phage)
  # rename target column, this allows using both clusters and KOs
  setnames(u.work, by, "target")
  u.mat <- dcast(u.work, phage ~ target, value.var = "u.phage")
  phages <- u.mat$phage
  u.mat[, phage := NULL]
  u.mat <- as.matrix(u.mat)
  rownames(u.mat) <- as.character(phages)
  X <- tcrossprod(u.mat)
  # eigenvals/vecs of the phage matrix
  phage.eig <- eigen(X)
  phage.centrality <- data.table(
    phage = rownames(X),
    centrality = abs(phage.eig$vectors[, 1])
  )
  prot.centrality <- t(u.mat) %*% phage.centrality$centrality
  prot.centrality <- data.table(target = rownames(prot.centrality),
                                centrality = prot.centrality[, 1])
  # restore column name to what was input
  setnames(prot.centrality, "target", by)
  return(list(
    phage = phage.centrality,
    protein = prot.centrality))
}

#' Creates randomized infection submatrix by permuting columns and/or rows and
#' calculates the generalized host and host protein JSDs
PermRowColNullModel <- function(ph, im, hm, cols = TRUE, rows = FALSE) {
  # get the infection submatrix for phages
  sm <- im[, ph]
  # permute rows and columns, preserving only total # infections
  if (cols) {
    # permute columns
    # this changes 'which hosts are infected'
    sm <- apply(sm, 2, function(v) {sample(v)})
  }
  if (rows) {
    # DIVERGENCE CALCULATIONS FAIL IF ANY COLS ARE ALL 0
    # so repeat until no cols sum to 0 (probably slow...sigh)
    test <- rep(0, ncol(sm))
    while (any(test == 0)) {
      # permuting rows requires transpose...silly
      # this changes 'which phages infect them'
      sm <- t(apply(sm, 1, function(v) {sample(v)}))
      test <- colSums(sm)
    }
  }
  # calculate host jsd
  host.jsd <- genJSD(sm)
  # make the new hostprot-copies-in-host range matrix
  ppmat <- hm %*% sm
  hostprot.jsd <- genJSD(ppmat)
  return(list(host.jsd, hostprot.jsd))
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