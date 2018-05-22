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

jsd <- fread("jsds/cholera.txt") # run make-jsds-cholera.R to precompute
js.dist <- dcast(jsd, sample.x ~ sample.y, value.var = "jsd")
rn <- js.dist$sample.x
js.dist <- sqrt(as.matrix(js.dist[, -1]))
rownames(js.dist) <- rn
js.dist <- js.dist[sample.names, sample.names]

# tdamapper method --------------------------------------------------------

source("k.first.R")
mds2 <- cmdscale(js.dist, 2)
po <- 70
ni <- c(10, 10)
mpr <- mapper2D(js.dist, filter_values = list(mds2[,1], mds2[,2]),
                percent_overlap = po,
                num_intervals = ni)

# characterize mapper vertices

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

# plot density
gordon.samples[, knn := kNN[sample]]
setkey(gordon.samples, diagnosis)
ggplot(gordon.samples["diarrhea"], aes(x = hour, y = knn)) +
  geom_point() +
  facet_wrap(~ subject, ncol = 2) +
  labs(y = "kNN", x = "hour", title = "diarrhea samples")

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
lo <- kk.fr(g1)
# lo <- create_layout(g1, layout = "igraph", algorithm = "fr", niter = 1000)
source("plot.mapper.R")
# color by fraction samples marked diarrhea
plot.mapper(lo, aes_(size = ~size, color = ~fd),
            list(size = "# samples", color = "f diarrhea")) +
  scale_color_distiller(palette = "Spectral")
# color vertices by mean knn density
plot.mapper(lo, aes_(size = ~size, color = ~mean.knn),
            list(size = "# samples", color = "mean kNN")) +
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

#' # 'Meta' persistent homology on the Mapper graph
#' ## b0 vs density threshold
# b0 vs density threshold ----------------------------------------------

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

#' ## Video of graph expansion as vertices are added from least to most dense
# global density level set video frames -----------------------------------
lvlset.grafs <- lapply(ls, function(th, graf, layout) {
  vids <- which(V(graf)$mean.knn <= th)
  sg <- induced_subgraph(graf, vids)
  xl <- c(min(V(graf)$x) - 1, max(V(graf)$x) + 1)
  yl <- c(min(V(graf)$y) - 1, max(V(graf)$y) + 1)
  cl <- sapply(c("min", "max"), function(fn, x) {
    do.call(fn, list(x= x))
    }, x = V(graf)$mean.knn)
  ggraph(sg, "manual", node.positions = data.frame(x = V(sg)$x, y = V(sg)$y)) +
    geom_edge_link0() +
    geom_node_point(aes(color = mean.knn, size = size)) +
    coord_cartesian(xlim = xl, ylim = yl) +
    scale_color_gradient2(midpoint = median(V(graf)$mean.knn),
                          high = "blue", mid = "yellow", low = "red",
                          limits = cl) +
    theme_graph()
}, graf = g1, layout = lo)
for (i in seq_along(ls)) {
  fn <- paste0("level-set/frames/", i, ".png")
  save_plot(fn, lvlset.grafs[[i]])
}

