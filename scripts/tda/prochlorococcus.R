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
jsds <- fread("jsds/prochlorococcus.txt")
jsds[, distance := sqrt(jsd)]
dist.mat <- dcast(jsds, sample.x ~ sample.y, value.var = "distance")
rn <- dist.mat[, sample.x]
dist.mat <- as.matrix(dist.mat[, -1])
rownames(dist.mat) <- rn
samples <- unique(prochlorococcus[, -c("ecotype", "abundance")])
samples[, phase := sin(pi * month / 2 / 12)]

# mds ---------------------------------------------------------------------

mds <- cmdscale(dist.mat, eig = TRUE)
mds$GOF # not too bad actually
plot(mds$points)

# mapper ------------------------------------------------------------------

po <- 50
ni <- c(20, 20)
mpr <- mapper2D(dist.mat, list(mds$points[, 1], mds$points[, 2]),
                percent_overlap = po, num_intervals = ni)
v2p <- vertex.2.points(mpr$points_in_vertex)
setnames(v2p, "point.name", "sample")
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
                    median.depth = median(depth, na.rm = TRUE),
                    median.temp = median(temp, na.rm = TRUE),
                    median.sal = median(sal, na.rm = TRUE),
                    median.phase = median(phase, na.rm = TRUE)
                    ),
                by = .(vertex, vertex.name)]
graf <- mapper.2.igraph(mpr)
# graf <- induced_subgraph(graf, V(graf)$size > 1)
setkey(vertices, vertex)
for (att in c("f.bats", "median.temp", "median.sal", "median.depth",
              "mean.phase")) {
  graf <- set_vertex_attr(graf, att, V(graf),
                          vertices[as.numeric(V(graf)), get(att)])
}
set.seed(1)
lo <- create_layout(graf, "fr", niter = 500)
plot.mapper(lo, aes_(size = ~size, color = ~f.bats)) +
  scale_color_distiller(palette = "Spectral")
plot.mapper(lo, aes_(size = ~size, color = ~median.temp)) +
  scale_color_distiller(palette = "Spectral")
plot.mapper(lo, aes_(size = ~size, color = ~median.sal)) +
  scale_color_distiller(palette = "Spectral")
plot.mapper(lo, aes_(size = ~size, color = ~median.depth)) +
  scale_color_distiller(palette = "Spectral")
