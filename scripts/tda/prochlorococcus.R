library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(cowplot)
library(combinat)

scripts.dir <- "../r/"
source(paste0(scripts.dir, "load-prochlorococcus-data.R"))
source("k.first.R")
source("vertex.2.points.R")
source("mapper.2.igraph.R")
source("plot.mapper.R")
source("k.first.R")
source("month.2.phase.R")
jsds <- fread("jsds/prochlorococcus.txt")
jsds[, distance := sqrt(jsd)]
dist.mat <- dcast(jsds, sample.x ~ sample.y, value.var = "distance")
rn <- dist.mat[, sample.x]
dist.mat <- as.matrix(dist.mat[, -1])
rownames(dist.mat) <- rn
colnames(dist.mat) <- rn
samples <- unique(prochlorococcus[, -c("ecotype", "abundance")])
samples[, phase := month.2.phase(cal.month)]

#' kNN density distribution is very skewed, but can smooth by increasing k:

# kNN ---------------------------------------------------------------------

k <- floor(nrow(samples) / 10)
knn <- apply(dist.mat, 1, k.first, k = k)
hist(knn)
set(samples, NULL, "knn", knn[samples$sample])

# l-infinity --------------------------------------------------------------

# linf <- linf(dist.mat)
# hist(linf)

# mds ---------------------------------------------------------------------

mds <- cmdscale(dist.mat, eig = TRUE) # 2D
#' GOF is good, but we find that samples are very unevenly distributed across
#' the 2D MDS-space, which is bad for Mapper:
mds$GOF # not too bad actually
plot(mds$points)
hist(mds$points[, 1])
hist(mds$points[, 2])
#' Converting to rank alleviates the problem somewhat.
#' Marginal distributions will be uniform, by definition.
rk.mds <- apply(mds$points, 2, rank, ties.method = "first")
plot(rk.mds)

# mapper ------------------------------------------------------------------

po <- 60
ni <- c(20, 20)
nb <- 10
ftr <- list(rk.mds[, 1], rk.mds[, 2])
mpr <- mapper2D(dist.mat, ftr,
                percent_overlap = po, num_intervals = ni,
                num_bins_when_clustering = nb)
v2p <- vertex.2.points(mpr$points_in_vertex)
v2p$sample <- rownames(dist.mat)[v2p$point]
setkey(v2p, sample)
setkey(samples, sample)
v2p <- samples[v2p]
v2p[, depth := as.numeric(depth)]
vertices <- v2p[, .(size = .N,
                    f.bats = sum(site == "bats") / .N,
                    mean.depth = mean(depth, na.rm = TRUE),
                    mean.temp = mean(temp, na.rm = TRUE),
                    mean.sal = mean(sal, na.rm = TRUE),
                    mean.phase = mean(phase, na.rm = TRUE),
                    mean.phase.bats = mean(phase[site == "bats"], na.rm = TRUE),
                    mean.phase.hot = mean(phase[site == "hot"], na.rm = TRUE),
                    mean.calmonth = mean(cal.month, na.rm = TRUE),
                    median.depth = median(depth, na.rm = TRUE),
                    median.temp = median(temp, na.rm = TRUE),
                    median.sal = median(sal, na.rm = TRUE),
                    mean.knn = mean(knn, na.rm = TRUE)
                    ),
                by = .(vertex, vertex.name)]
graf <- mapper.2.igraph(mpr)
# graf <- induced_subgraph(graf, V(graf)$size > 1)
setkey(vertices, vertex)
for (att in names(vertices)[-c(1, 2)]) {
  graf <- set_vertex_attr(graf, att, V(graf),
                          vertices[as.numeric(V(graf)), get(att)])
}
set.seed(1)
lo <- create_layout(graf, "fr", niter = 500)
#' # Composition varies continuously with temperature
plot.mapper(lo, aes_(size = ~size, color = ~mean.temp)) +
  scale_color_distiller(palette = "Spectral")
#' # Composition varies continuously with depth
plot.mapper(lo, aes_(size = ~size, color = ~mean.depth)) +
  scale_color_distiller(palette = "Spectral")
#' # Composition is more stable at low depth and high temperature
plot.mapper(lo, aes_(size = ~size, color = ~mean.knn)) +
  scale_color_distiller(palette = "Spectral")
#' # Composition is not well-separated by site
plot.mapper(lo, aes_(size = ~size, color = ~f.bats)) +
  scale_color_distiller(palette = "Spectral")
#' # Composition is not well-separated by salinity
#'
#' BATS is also systematically higher salinity than HOT, so this is skewed.
plot.mapper(lo, aes_(size = ~size, color = ~mean.sal)) +
  scale_color_distiller(palette = "Spectral")
#' # Composition is not well-separated by month
#' Letting March be phase = 0.
plot.mapper(lo, aes_(size = ~size, color = ~mean.phase)) +
  scale_color_distiller(palette = "Spectral")
plot.mapper(lo, aes_(size = ~size, color = ~mean.phase.bats)) +
  scale_color_distiller(palette = "Spectral")
plot.mapper(lo, aes_(size = ~size, color = ~mean.phase.hot)) +
  scale_color_distiller(palette = "Spectral")