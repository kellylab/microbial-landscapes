library(data.table)
library(TDAmapper)
library(philentropy)
library(cowplot)

util.dir <- "../r/"

source(paste0(util.dir, "load_cholera_data.R"))

# load data ---------------------------------------------------------------

gordon[, freq := count / sum(count), by = sample]
distribs <- dcast(gordon, sample ~ otu, value.var = "freq", fill = 0)
sample.names <- distribs$sample
distribs <- as.matrix(distribs[, -1])

# compute js distance -----------------------------------------------------

jsd <- JSD(distribs)
rownames(jsd) <- sample.names
colnames(jsd) <- sample.names
js.dist <- sqrt(jsd)

# tdamapper method --------------------------------------------------------

# filter: mean distance to k nearest pts
k.nearest <- function(v, k) {
  v <- sort(v[v > 0])
  mean(v[1:k])
}
filter.position.knn <- function(dst, k) {
  f1 <- apply(dst, 1, k.nearest, k = k) # filtering by local density
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
mpr <- mapper2D(js.dist, filter_values = ftr, percent_overlap = 50,
                num_intervals = c(5, 5))

# characterize mapper vertices --------------------------------------------


library(igraph)
library(ggraph)
g1 <- graph.adjacency(mpr$adjacency, mode="undirected")
# size
V(g1)$size <- sapply(mpr$points_in_vertex, length)
# fraction diarrhea
fstate <- function(v, dt) {
  setkey(dt, sample)
  snames <- names(v)
  s <- dt[snames, state]
  sum(s == "diarrhea") / length(s)
}
V(g1)$fd <- sapply(mpr$points_in_vertex, fstate, dt = gordon)
# density
vert.knn <- function(ps, js.dist, k, agg = "mean") {
  knn <- sapply(ps, function(p, js.dist, k) {
    v <- js.dist[p,]
    k.nearest(v, k)
  }, js.dist = js.dist, k = k)
  do.call(agg, list(x = knn))
}
V(g1)$mean.knn <- sapply(mpr$points_in_vertex, vert.knn, js.dist = js.dist,
                         k = 10)
# time
V(g1)$mean.t <- sapply(mpr$points_in_vertex, function(p, dt) {
  setkey(dt, sample)
  snames <- names(p)
  dt[snames, mean(hour)]
}, dt = gordon)


# draw graphs -------------------------------------------------------------


kk.fr <- function(graf) {
  l1 <- create_layout(graf, layout = "igraph", algorithm = "kk")
  l1 <- create_layout(graf, layout = "igraph", algorithm = "fr",
                      niter = 1000,
                      coords = as.matrix(l1[, c("x", "y")]))
  l1
}
lo <- kk.fr(g1)
ggraph(lo) +
  geom_edge_link() +
  geom_node_point(aes(size = size, color = fd)) +
  labs(size = "# samples", color = "f diarrhea") +
  theme(aspect.ratio = 1) +
  scale_color_distiller(palette = "Spectral") +
  theme_graph(base_family = "Helvetica")
# color vertices by mean knn density
ggraph(lo) +
  geom_edge_link() +
  geom_node_point(aes(size = size, color = mean.knn)) +
  labs(size = "# samples") +
  theme(aspect.ratio = 1) +
  scale_color_distiller(palette = "Blues") +
  theme_graph(base_family = "Helvetica")

ggraph(lo) +
  geom_edge_link() +
  geom_node_point(aes(size = size, color = mean.t)) +
  labs(size = "# samples") +
  theme(aspect.ratio = 1) +
  scale_color_distiller(palette = "Spectral") +
  theme_graph(base_family = "Helvetica")

# subject trajectories ----------------------------------------------------
vtxmap <- lapply(mpr$points_in_vertex, function(v) {
  data.table(sample = names(v))
})
vtxmap <- rbindlist(vtxmap, idcol = "vertex")
V(g1)$x <- lo[V(g1), "x"]
V(g1)$y <- lo[V(g1), "y"]
sample.spx <- lapply(sample.names, function(s, vtxmap, g) {
  setkey(vtxmap, sample)
  vs <- vtxmap[s, vertex]
  if (length(vs) == 1) {
    if (is.na(vs)) {
      NA
    } else {
      induced_subgraph(g, vs)
    }
  } else{
    induced_subgraph(g, vs)
  }
}, vtxmap = vtxmap, g = g1)
names(sample.spx) <- sample.names
sample.spx <- sample.spx[!is.na(sample.spx)]
graph.2.df <- function(g) {
  data.frame(x = V(g)$x, y = V(g)$y, size = V(g)$size)
}
sample.overlays <- lapply(names(sample.spx), function(name, ss, lo) {
  sg <- sample.spx[[name]]
  if (grepl("diarrhea", name)) {
    clr <- "red"
  } else {
    clr <- "blue"
  }
  df <- graph.2.df(sg)
  ggraph(lo) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey50") +
    theme_graph(base_family = "Helvetica") +
    labs(title = name) +
    geom_point(aes(x = x, y = y, size = size),
               data = graph.2.df(sg),
               color = clr)
}, ss = sample.spx, lo = lo)
names(sample.overlays) <- names(sample.spx)
for (s in names(sample.overlays)) {
  ggsave(paste0("subject-trajectories/", s, ".png"),
         sample.overlays[[s]], scale = 0.5)
}


# ‘meta’ persistent homology ----------------------------------------------

V(g1)$max.knn <- sapply(mpr$points_in_vertex, vert.knn, js.dist = js.dist,
                        k = 10, agg = "max")
ls <- seq(min(V(g1)$max.knn), max(V(g1)$max.knn), length.out = 1000)
ls.gs <- lapply(ls, function(x, g) induced_subgraph(g, V(g)$max.knn <= x),
                g = g1)
ncomps <- sapply(ls.gs, function(g) components(g)$no)
ggplot(data.frame(pi = ls, b0 = ncomps), aes(x = pi, y = b0)) +
  geom_line()
rep.vs <- which(V(g1)$max.knn <= max(ls[ncomps == max(ncomps)]))
# rep.vs <- which(V(g1)$mean.knn <= ls[60])
rep.xy <- lo[rep.vs, c("x", "y")]
rep.g <- induced_subgraph(g1, rep.vs)
ggraph(rep.g, layout = "manual", node.positions = rep.xy) +
  geom_edge_link() +
  geom_node_point(aes(size = size, color = max.knn)) +
  labs(size = "# samples") +
  theme(aspect.ratio = 1) +
  scale_color_distiller(palette = "Blues") +
  theme_graph(base_family = "Helvetica")
ggraph(rep.g, layout = "manual", node.positions = rep.xy) +
  geom_edge_link() +
  geom_node_point(aes(size = size, color = fd)) +
  labs(size = "# samples") +
  theme(aspect.ratio = 1) +
  scale_color_distiller(palette = "Spectral") +
  theme_graph(base_family = "Helvetica")
