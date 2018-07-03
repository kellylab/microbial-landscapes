library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(tidygraph)
library(cowplot)

scripts.dir <- "../r/"
source(paste0(scripts.dir, "load-nahant-data.R"))
source("k.first.R")
source("vertex.2.points.R")
source("mapper.2.igraph.R")
source("plot.mapper.R")

nahant <- nahant[!is.na(kingdom)]
nahant[, freq := value / sum(value), by = .(kingdom, day)]

# jsds

jsds <- lapply(split(nahant, by = "kingdom"), function(dt) {
  x <- dcast(dt, day ~ OTU, value.var = "freq", fill = 0)
  rn <- x$day
  x <- as.matrix(x[, -1])
  rownames(x) <- rn
  m <- JSD(x)
  rownames(m) <- rn
  colnames(m) <- rn
  m
})
distances <- lapply(jsds, sqrt)
k <- 10
samples <- lapply(distances, function(m, k) {
  v <- apply(m, 1, k.first, k = k)
  data.table(day = as.numeric(names(v)), kNN = v)
}, k = k)
#' Bacteria show trend of increasing instability toward day 250, while
#' eukaryotes have about the same stability:
ggplot(rbindlist(samples, idcol = "kingdom"), aes(x = day, y = kNN)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ kingdom, nrow = 2)
mds <- lapply(distances, cmdscale, eig = TRUE)
#' 2D MDS is lossy, but suggests existence of 2 bacterial and 2-3 eukaryotic
#' states:
for (x in mds) {
  print(x$GOF)
  plot(x$points)
}
rk.mds <- lapply(mds, function(x) {
  pts <- x$points
  apply(pts, 2, rank, ties.method = "first")
})


# mapper ------------------------------------------------------------------

ni <- c(10, 10)
po <- 70
mpr <- mapply(function(distance, rk.mds) {
  mapper2D(distance, list(rk.mds[, 1], rk.mds[, 2]), num_intervals = ni,
           percent_overlap = po)
}, distance = distances, rk.mds = rk.mds, SIMPLIFY = FALSE)
v2p <- lapply(mpr, function(x) {
  dt <- vertex.2.points(x$points_in_vertex)
  setnames(dt, "point.name", "day")
  dt[, day := as.numeric(day)]
  dt
})
v2p <- mapply(merge, x = v2p, y = samples, MoreArgs = list(by = "day"),
              SIMPLIFY = FALSE)
vertices <- lapply(v2p, function(dt) {
  dt[, .(mean.day = mean(day),
         mean.kNN = mean(kNN),
         size = .N),
     by = .(vertex, vertex.name)]
})
grafs <- lapply(mpr, mapper.2.igraph)
grafs <- mapply(function(graf, verts) {
  graf %>%
    as_tbl_graph %>%
    activate(nodes) %>%
    left_join(verts, by = c("name" = "vertex.name"))
}, graf = grafs, verts = vertices, SIMPLIFY = FALSE)
set.seed(0)
layouts <- lapply(grafs, create_layout, layout = "fr")
size.plots <- mapply(function(lo, kingdom) {
  p <- plot.mapper(lo, aes_(size = ~size, color = ~mean.day),
                   list(title = kingdom))
  p + scale_color_distiller(palette = "Spectral") +
    theme(legend.position = "bottom")
}, lo = layouts, kingdom = names(layouts), SIMPLIFY = FALSE)
plot_grid(plotlist = size.plots)
kNN.plots <- mapply(function(lo, kingdom) {
  p <- plot.mapper(lo, aes_(size = ~size, color = ~mean.kNN),
                   list(title = kingdom))
  p + scale_color_distiller(palette = "Spectral") +
    theme(legend.position = "bottom")
}, lo = layouts, kingdom = names(layouts), SIMPLIFY = FALSE)
plot_grid(plotlist = kNN.plots)


# frames ------------------------------------------------------------------

grafs <- mapply(function(graf, lo) {
  graf %>%
    activate(nodes) %>%
    mutate(x = lo$x, y = lo$y)
}, graf = grafs, lo = layouts, SIMPLIFY = FALSE)
sample.subgraphs <- mapply(function(graf, v2p) {
  spl <- split(v2p, by = "day")
  lapply(spl, function(s, graf) {
    graf %>%
      slice(s$vertex) %>%
      mutate(day = s$day)
  }, graf = graf)
}, graf = grafs, v2p = v2p, SIMPLIFY = FALSE)
sample.plots <- mapply(function(sgs, layout) {
  lapply(sgs, function(sg, layout) {
    sg <- as.data.frame(sg)
    day <- unique(sg$day)
    p <- plot.mapper(layout, aes_(size = ~size), NULL, color = "grey")
    p + geom_point(aes(x = x, y = y, size = size),
                 data = sg, color = "blue") +
      labs(title = day)
  }, layout = layout)
}, sgs = sample.subgraphs, layout = layouts)
for (dt in samples) {
  setorder(dt, day)
  dt[, id := seq_len(.N)]
}
for (kingdom in names(samples)) {
  for (x in seq_len(nrow(samples[[1]]))) {
    day <- samples[[kingdom]][x, day]
    fn <- paste0(kingdom, x, ".png")
    save_plot(paste0("nahant-trajectories/frames/", fn),
              sample.plots[[kingdom]][[x]])
  }
}
