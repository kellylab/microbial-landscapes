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

# density estimator on a grid method
# mds
#n <- length(rn)
#k <- 2
#x <- cmdscale(js.dist, k = k, eig = TRUE)
#print("MDS goodness of fit:")
#print(x$GOF)
#mds <- x$points
#ngrid <- 100
#grid <- expand.grid(
#  as.data.frame(
#  apply(mds, 2, function(v) seq(min(v), max(v), length.out = ngrid)
#)))

# create filtration -------------------------------------------------------

#maxscale <- max(js.dist)
#print(paste("Maxscale = max distance =" , maxscale))
#
#rips <- ripsFiltration(js.dist, maxdimension = 0, maxscale = maxscale,
#  dist = "arbitrary", library = "Dionysus")

# pdf("figs/cholera_persistence_diag.pdf")
# plot(rips$diagram, main = "Diagram")
# plot(rips$diagram, barcode = TRUE, main = "Barcode")
# cdf <- ecdf(rips$diagram[, 3])
# plot(cdf, main = "b0")
# bzero <- function(x) (1 - cdf(x)) * nrow(rips$diagram)
# plot(bzero, main = "# components", log = "y")
# dev.off()

# tdamapper method --------------------------------------------------------

# filter: mean distance to k nearest pts
k.nearest <- function(v, k) {
  v <- sort(v[v > 0])
  mean(v[1:k])
}
f1 <- apply(js.dist, 1, k.nearest, k = 10)
setkey(gordon, sample)
f2 <- gordon[sample.names, unique(hour)]
mpr <- mapper1D(js.dist, filter_values = f1, num_intervals = 100, percent_overlap = 90)
#mpr <- mapper2D(js.dist, filter_values = list(f1, f2))
library(igraph)
library(ggraph)
g1 <- graph.adjacency(mpr$adjacency, mode="undirected")
V(g1)$size <- sapply(mpr$points_in_vertex, length)
fstate <- function(v, dt) {
  setkey(dt, sample)
  snames <- names(v)
  s <- dt[snames, state]
  sum(s == "diarrhea") / length(s)
}
V(g1)$fd <- sapply(mpr$points_in_vertex, fstate, dt = gordon)
#plot(g1, layout=layout.auto(g1))
kk.fr <- function(graf) {
  l1 <- create_layout(graf, layout = "igraph", algorithm = "kk")
  l2 <- create_layout(graf, layout = "igraph", algorithm = "fr",
                      coords = as.matrix(l1[, c("x", "y")]))
  l2
}
lo <- kk.fr(g1)
ggraph(lo) +
  geom_edge_link() +
  geom_node_point(aes(size = size, color = fd)) +
  labs(size = "# samples", color = "f diarrhea") +
  theme(aspect.ratio = 1) +
  scale_color_distiller(palette = "Spectral") +
  # scale_size_continuous(range = c(0.1, 1)) +
  theme_graph(base_family = "Helvetica")