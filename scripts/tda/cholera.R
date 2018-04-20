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
# ftr <- filter.position2D(js.dist)
ftr <- filter.position.knn(js.dist, 10)
po <- 50
ni <- c(5, 5)
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
lo <- kk.fr(g1)
source("plot.mapper.R")
# color by fraction samples marked diarrhea
plot.mapper(lo, aes_(size = ~size, color = ~fd),
            list(size = "# samples", color = "f diarrhea")) +
  scale_color_distiller(palette = "Spectral")
# color vertices by mean knn density
plot.mapper(lo, aes_(size = ~size, color = ~mean.knn),
            list(size = "# samples")) +
  scale_color_distiller(palette = "Blues")
# color by mean time in estimated hours
plot.mapper(lo, aes_(size = ~size, color = ~mean.t), list(size = "samples")) +
  scale_color_distiller(palette = "Spectral")

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

V(g1)$max.knn <- sapply(mpr$points_in_vertex, vertex.attribute,
                        point.attributes = kNN, summ = "max")
ls <- sort(unique(V(g1)$max.knn))
ls.gs <- lapply(ls, function(x, g) induced_subgraph(g, V(g)$max.knn <= x),
                g = g1)
comps <- lapply(ls.gs, components)
ncomps <- sapply(comps, function(g) g$no)
ggplot(data.frame(pi = ls, b0 = ncomps), aes(x = pi, y = b0)) +
  geom_line() + geom_point()

#' Getting the representative level set as the one with most components.
#' Get the vertices in that level set:
sel <- max(which(ncomps == max(ncomps))) # level set index
rep.vs <- which(V(g1)$max.knn <= ls[sel])
rep.xy <- lo[rep.vs, c("x", "y")]
rep.g <- induced_subgraph(g1, rep.vs)
V(rep.g)$state <- as.character(comps[[sel]]$membership)
rep.lo <- create_layout(rep.g, "manual", node.positions = rep.xy)
plot.mapper(rep.lo, aes_(size = ~size, color = ~max.knn),
            list(size = "# samples")) +
  scale_color_distiller(palette = "Blues")
plot.mapper(rep.lo, aes_(size = ~size, color = ~fd),
            list(size = "# samples")) +
  scale_color_distiller(palette = "Spectral")
plot.mapper(rep.lo, aes_(size = ~size, color = ~state),
            list(size = "# samples"))

#' Characterize each of the components in the representative level set.
state.knns <- sapply(unique(V(rep.g)$state), function(si, graph, fn = "max") {
  sg <- induced_subgraph(graph, V(graph)$state == si)
  knn <- V(sg)$max.knn
  c(max = max(knn), min = min(knn))
}, graph = rep.g)
state.knns <- as.data.table(t(state.knns))
state.knns$state <- factor(seq(nrow(state.knns)),
                           levels = as.character(seq(nrow(state.knns))))
ggplot(state.knns, aes(x = min, y = state)) +
  geom_segment(aes(group = state, xend = max, yend = state)) +
  geom_point(data = function(dt) dplyr::filter(dt, max == min)) +
  labs(x = "max kNN", color = "persistence")

#' Subject time points characterized by state(s)
state.knns$persistence <- state.knns$max - state.knns$min
state.knns[, rank := frank(-persistence)]

#' Get the vertices associated with the top 3 states
stable.vertices <- lapply(state.knns[rank <= 3, state], function(st, rg) {
  data.table(state = st, vertex = V(rg)$name[V(rg)$state == st])
}, rg = rep.g)
stable.vertices <- rbindlist(stable.vertices)
setkey(gordon.samples, sample)
setkey(vtxmap, sample)
vtxmap <- vtxmap[gordon.samples]
setkey(stable.vertices, vertex)
setkey(vtxmap, vertex.name)
vtxmap <- stable.vertices[vtxmap]
ggplot(vtxmap, aes(x = idx, y = state)) +
  geom_point(aes(color = diagnosis)) +
  theme(panel.grid.major.y = element_line(size = 0.1, color = "grey50")) +
  facet_wrap(~ subject, scales = "free_x")

#' Just the diarrhea samples, since the time scale changes in recovery
setkey(vtxmap, diagnosis)
ggplot(vtxmap["diarrhea"], aes(x = hour, y = state)) +
  geom_point() +
  theme(panel.grid.major.y = element_line(size = 0.1, color = "grey50")) +
  labs(title = "diarrhea") +
  facet_wrap(~ subject)
