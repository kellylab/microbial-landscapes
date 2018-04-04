library(data.table)
library(TDA)
library(philentropy)

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
vert.knn <- function(ps, js.dist, k) {
  knn <- sapply(ps, function(p, js.dist, k) {
    v <- js.dist[p,]
    k.nearest(v, k)
  }, js.dist = js.dist, k = k)
  mean(knn)
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
