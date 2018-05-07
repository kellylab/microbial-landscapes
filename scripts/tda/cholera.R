library(data.table)
library(TDAmapper)
library(philentropy)
library(cowplot)
library(tidyr)

util.dir <- "../r/"

source(paste0(util.dir, "load_cholera_data.R"))

# load data ---------------------------------------------------------------

gordon[, freq := count / sum(count), by = sample]
gordon.samples <- unique(gordon[, .(sample, subject, diagnosis, id, hour)])
gordon.samples[, idx := frank(hour), by = subject]
distribs <- dcast(gordon, sample ~ otu, value.var = "freq", fill = 0)
sample.names <- distribs$sample
distribs <- as.matrix(distribs[, -1])

# compute js distance -----------------------------------------------------

jsd <- JSD(distribs)
rownames(jsd) <- sample.names
colnames(jsd) <- sample.names
js.dist <- sqrt(jsd)

# tdamapper method --------------------------------------------------------

source("k.first.R")
filter.position.knn <- function(dst, k) {
  f1 <- apply(dst, 1, k.first, k = k) # filtering by local density
  mds1 <- cmdscale(js.dist, 1) # filtering by 'position'
  f2 <- mds1[, 1]
  list(f1, f2)
}
filter.position2D <- function(dst) {
  # filtering by position only
  mds2 <- cmdscale(dst, 2)
  list(mds2[, 1], mds2[, 2])
}
ftr <- filter.position2D(js.dist)
# ftr <- filter.position.knn(js.dist, 10)
po <- 80
ni <- c(10, 10)
mpr <- mapper2D(js.dist, filter_values = ftr, percent_overlap = po,
                num_intervals = ni)
# fix the adjacency matrix
source("mapper.adj.R")
mpr$adjacency <- mapper.adj(mpr$points_in_vertex)

# characterize mapper vertices --------------------------------------------

library(igraph)
library(ggraph)
g1 <- graph.adjacency(mpr$adjacency, mode="undirected")
# name vertices so can track to subgraphs later
V(g1)$name <- paste0("v", seq_along(V(g1)))
# size
V(g1)$size <- sapply(mpr$points_in_vertex, length)
# fraction diarrhea
source("named.vector.R")
sample.diagnosis <- named.vector(gordon.samples$diagnosis,
                                 gordon.samples$sample)
fstate <- function(x) {
  sum(x == "diarrhea") / length(x)
}
source("vertex.attribute.R")
V(g1)$fd <- sapply(mpr$points_in_vertex, vertex.attribute,
                   point.attributes = sample.diagnosis, summ = "fstate")
# density
k <- 10
kNN <- apply(js.dist, 1, k.first, k = k)
V(g1)$mean.knn <- sapply(mpr$points_in_vertex, vertex.attribute,
                         point.attributes = kNN)
# time
sample.times <- named.vector(gordon.samples$hour, gordon.samples$sample)
V(g1)$mean.t <- sapply(mpr$points_in_vertex, vertex.attribute,
                       point.attributes = sample.times)

# draw graphs -------------------------------------------------------------

set.seed(1)
kk.fr <- function(graf) {
  l1 <- create_layout(graf, layout = "igraph", algorithm = "kk")
  l1 <- create_layout(graf, layout = "igraph", algorithm = "fr",
                      niter = 1000,
                      coords = as.matrix(l1[, c("x", "y")]))
  l1
}
# lo <- kk.fr(g1)
lo <- create_layout(g1, layout = "igraph", algorithm = "fr", niter = 1000)
source("plot.mapper.R")
# color by fraction samples marked diarrhea
plot.mapper(lo, aes_(size = ~size, color = ~fd),
            list(size = "# samples", color = "f diarrhea")) +
  scale_color_distiller(palette = "Spectral")
# color vertices by mean knn density
plot.mapper(lo, aes_(size = ~size, color = ~mean.knn),
            list(size = "# samples")) +
  scale_color_distiller(palette = "Spectral", direction = 1)
# color by mean time in estimated hours
plot.mapper(lo, aes_(size = ~size, color = ~mean.t), list(size = "samples")) +
  scale_color_distiller(palette = "Spectral", direction = 1)

# subject trajectories ----------------------------------------------------
source("vertex.2.points.R")
vtxmap <- vertex.2.points(mpr$points_in_vertex)
setnames(vtxmap, "point.name", "sample")
V(g1)$x <- lo[V(g1), "x"]
V(g1)$y <- lo[V(g1), "y"]
# produce a subgraph for every sample
source("sample.subgraphs.R")
sample.spx <- sample.subgraphs(
  vtxmap[, .(samples = sample, vertices = vertex)], g1)
sample.spx <- sample.spx[!is.na(sample.spx)]
sample.overlays <- lapply(names(sample.spx), function(name, ss, lo) {
  sg <- sample.spx[[name]]
  if (grepl("diarrhea", name)) {
    clr <- "red"
  } else {
    clr <- "blue"
  }
  plot.mapper(lo, aes_(size = ~size), list(title = name), color = "grey50") +
    geom_point(aes(x = x, y = y, size = size),
               data = as.data.frame(vertex_attr(sg)), color = clr)
}, ss = sample.spx, lo = lo)
names(sample.overlays) <- names(sample.spx)
for (s in names(sample.overlays)) {
  save_plot(paste0("subject-trajectories/frames/", s, ".png"),
            sample.overlays[[s]])
}

# ‘meta’ persistent homology ----------------------------------------------

V(g1)$mean.knn <- sapply(mpr$points_in_vertex, vertex.attribute,
                        point.attributes = kNN, summ = "mean")
ls <- sort(unique(V(g1)$mean.knn))
ls.gs <- lapply(ls, function(x, g) induced_subgraph(g, V(g)$mean.knn <= x),
                g = g1)
comps <- lapply(ls.gs, components)
# ncomps <- sapply(comps, function(g) sum(g$csize > 1)) # filter singletons
ncomps <- sapply(comps, function(g) g$no)
ggplot(data.frame(threshold = ls, b0 = ncomps), aes(x = threshold, y = b0)) +
  geom_line() + geom_point()

#' Get the longest interval in density threshold where b0 is constant.
#' First get the number of intervals in the threshold for each run where b0 is
#' constant:
b0.runs <- rle(ncomps)
b0.runs <- data.table(length = b0.runs$lengths, b0 = b0.runs$values)
#' Cumulatively sum over intervals in the threshold value to get the end
#' threshold value of each run:
b0.runs[, end := sapply(cumsum(length), function(n, th) {
  th[n]
}, th = ls)]
#' Shift the end values down and prepend 0 for the start threshold value of
#' each run:
b0.runs[, start := c(0, end[-.N])]
#' Get persistence length of each run in the dimension of the threshold.
b0.runs[, persistence := end - start]
#' Get rid of the run between 0 and the first point.
b0.runs <- b0.runs[-1,]
#' Find the threshold at the end of the run with maximum persistence
thc <- b0.runs[b0 > 1, end[persistence == max(persistence)]]

#' Get the vertices in that level set:
rep.vs <- which(V(g1)$mean.knn <= thc)
rep.xy <- lo[rep.vs, c("x", "y")]
rep.g <- induced_subgraph(g1, rep.vs)
V(rep.g)$state <- as.character(comps[[which(ls == thc)]]$membership)
rep.lo <- create_layout(rep.g, "manual", node.positions = rep.xy)
plot.mapper(rep.lo, aes_(size = ~size, color = ~state), list(size = "samples"))
# plot.mapper(lo, aes_(size = ~size), list(size = "samples"), color = "grey50") +
#   geom_point(data = rep.lo,
#              aes_(x = ~x, y = ~y, size = ~size, color = ~state))
