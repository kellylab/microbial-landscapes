library(philentropy)
library(TDAmapper)
library(data.table)
library(tidyverse)
library(ggraph)
library(igraph)
library(cowplot)

# setup -------------------------------------------------------------------


scripts.dir <- "../r/"
source(paste0(scripts.dir, "load_david_data.R"))
samples <- unique(david[, .(sample, subject, day)])

# if haven't pre-computed jsds, run make-jsds-david.R
jsds <- fread("jsds/david.txt")
jsds[, distance := sqrt(jsd)]
dist.mat <- dcast(jsds, sample.x ~ sample.y, value.var = "jsd")
sample.names <- dist.mat$sample.x
dist.mat <- as.matrix(dist.mat[, -1])
rownames(dist.mat) <- sample.names

#' # Preview of data using 2D MDS
#'
#' 2D MDS shows much more variance in dimension 1 than dimension 2, resulting in
#' large empty spaces.
#' This suggests that 2D MDS would not be the best choice for both filters.
mds2d <- cmdscale(dist.mat, 2)
plot(mds2d)

#' # Mapper
#' ## Filter by 1D MDS

# mds1 filter -------------------------------------------------------------


mds1d <- cmdscale(dist.mat, 1)
po <- 80
ni <- 50
mpr <- mapper1D(dist.mat, mds1d, percent_overlap = 70)
mpr$points_in_vertex <- lapply(mpr$points_in_vertex, function(v, dt) {
  names(v) <- dt[v, sample]
  v
}, dt = samples)

#' ## Filter by 2D MDS

# mds2 filter -------------------------------------------------------------



po <- 95
ni <- c(10, 2)
mpr <- mapper2D(dist.mat, list(mds2d[, 1], mds2d[, 2]), num_intervals = ni,
                percent_overlap = po)

# mds1 + kNN filter -------------------------------------------------------

mds1d <- cmdscale(dist.mat, 1)
source("k.first.R")
k <- 20
knn <- apply(dist.mat, 1, k.first, k = k)
po <- 80
ni <- c(5, 10)
mpr <- mapper2D(dist.mat, list(mds1d, knn), num_intervals = ni,
                percent_overlap = po)
mpr$points_in_vertex <- lapply(mpr$points_in_vertex, function(v, dt) {
  names(v) <- dt[v, sample]
  v
}, dt = samples)


# format mapper output ----------------------------------------------------


source("mapper.adj.R")
mpr$adjacency <- mapper.adj(mpr$points_in_vertex) # correct adj matrix

#' Format Mapper output:
source("vertex.2.points.R")
v2p <- vertex.2.points(mpr$points_in_vertex) # map vertices to samples
source("mapper.2.igraph.R")
graf <- mapper.2.igraph(mpr)
V(graf)$subject <- sapply(mpr$points_in_vertex, function(pts, dt) {
  pnames <- names(pts)
  setkey(dt, sample)
  dt[pnames, sum(subject == "A") / .N]
}, dt = samples)


#' ## Plots

# plot --------------------------------------------------------------------


source("plot.mapper.R")
#' Fraction of samples in each vertex belonging to each subject:
set.seed(2)
lo <- create_layout(graf, layout = "igraph", algorithm = "kk")
lo <- create_layout(graf, layout = "igraph", algorithm = "fr",
                    niter = 1000,
                    coords = as.matrix(lo[, c("x", "y")]))
plot.mapper(lo, aes_(size = ~size, color = ~subject),
            list(color = "fraction A")) +
  scale_color_gradient2(midpoint = 0.5, mid = "yellow")
