library(philentropy)
library(TDAmapper)
library(data.table)
library(tidyverse)
library(ggraph)
library(igraph)
library(tidygraph)
library(cowplot)

# setup -------------------------------------------------------------------


source("dist2knn.R")
source("vertex.2.points.R")
source("mapper.2.igraph.R")
source("plot.mapper.R")
scripts.dir <- "../r/"
figs.dir <- "../../figures/tda/"
source(paste0(scripts.dir, "load_david_data.R"))

samples <- unique(david[, .(sample, subject, day)])
samples[, day := as.numeric(day)]

# if haven't pre-computed jsds, run make-jsds-david.R
jsds <- fread("jsds/david.txt")
jsds[, distance := sqrt(jsd)]
dist.mat <- dcast(jsds, sample.x ~ sample.y, value.var = "jsd")
sample.names <- dist.mat$sample.x
dist.mat <- as.matrix(dist.mat[, -1])
rownames(dist.mat) <- sample.names

# calculate k nearest neighbors
k <- 10
kNN <- dist2knn(dist.mat, k)
samples[, kNN := kNN[sample]]

# event metadata
setkey(samples, subject)
events <- mapply(function(start, end, subject, event) {
  x <- data.table(day = seq(start, end, by = 1), subject = subject,
                  event = event)
  x
}, start = c(0, 71, 80, 104, 123, 0, 151, 160 ),
end = c(70, 122, 85, 113, samples["A", max(day)],
          150, 159, samples["B", max(day)]),
subject = c("A", "A", "A", "A", "A",
            "B", "B", "B"),
event = c("US (pre)", "travel", "diarrhea 1", "diarrhea 2", "US (post)",
          "pre-Salmonella", "Salmonella", "post-Salmonella"),
SIMPLIFY = FALSE) %>% rbindlist(use.names = TRUE)
# collapse event labels per day
events <- events[, .(event = paste(event, collapse = " + ")),
                 by = .(subject, day)]
samples <- merge(samples, events, by = c("subject", "day"))

#' # Data preview
#' ## Subjects traverse low density regions of state space near events
ggplot(samples, aes(x = day, y = kNN)) +
  geom_point(aes(color = event)) +
  stat_smooth(span = 0.25) +
  scale_y_log10() +
  facet_wrap(~ subject, ncol = 1)

#' ## MDS
#'
#' 2D MDS shows uneven distribution. Transform to rank to make more even.

# mds ------------------------------------------------------------------

mds2d <- cmdscale(dist.mat, eig = TRUE)
mds2d$GOF # not very good
plot(mds2d$points)
rk.mds <- apply(mds2d$points, 2, rank, ties.method = "first")
plot(rk.mds)

#' # Mapper
#' ## Filter by 2D MDS

# mds filter -------------------------------------------------------------

ftr <- list(rk.mds[, 1], rk.mds[, 2])
ni <- c(20, 20)
po <- 70
mpr <- mapper2D(dist.mat, ftr, num_intervals = ni, percent_overlap = po)
v2p <- vertex.2.points(mpr$points_in_vertex)
v2p <- merge(v2p, samples, by.x = "point.name", by.y = "sample")
vertices <- v2p[, .(
  subject = sum(subject == "A") / .N,
  mean.knn = mean(kNN),
  size = .N
), by = .(vertex, vertex.name)]

graf <- mapper.2.igraph(mpr)
graf <- graf %>%
  as_tbl_graph %>%
  activate(nodes) %>%
  left_join(vertices, by = c("name" = "vertex.name"))

#' ## Plots

#' Fraction of samples in each vertex belonging to each subject:
set.seed(0)
lo <- create_layout(graf, "fr")
plot.mapper(lo, aes_(size = ~size, color = ~subject),
            list(color = "fraction A")) +
  scale_color_distiller(palette = "Spectral")
save_plot(paste0(figs.dir, "david-fsubject.pdf"), last_plot(), base_height = 6)
plot.mapper(lo, aes_(size = ~size, color = ~mean.knn)) +
  scale_color_distiller(trans = "log10")


# subject frames ----------------------------------------------------------


graf <- graf %>%
  activate(nodes) %>%
  mutate(x = lo$x, y = lo$y)
sample.subgraphs <- lapply(split(v2p, by = "point.name"), function(s, graf) {
  graf %>%
    slice(s$vertex) %>%
    mutate(event = s$event, subject = s$subject, kNN = s$kNN,
           sample = s$point.name)
}, graf = graf)
sample.plots <- lapply(sample.subgraphs, function(sg, graf, events) {
  np <- as.data.frame(select(graf, x, y))
  sg <- as.data.frame(sg)
  sample <- unique(sg$sample)
  lo <- create_layout(graf, "manual", node.positions = np)
  p <- plot.mapper(lo, aes_(size = ~size), NULL, color = "grey")
  p + geom_point(aes(x = x, y = y, size = size, color = event),
                 data = sg) +
    scale_color_discrete(limits = events) +
    labs(title = sample)
}, graf = graf, events = unique(events$event))
setorder(samples, subject, day)
samples[, i := seq_len(.N), by = subject]
for (x in seq_len(nrow(samples))) {
  s <- samples[x, sample]
  sj <- samples[x, subject]
  fn <- paste0(sj, samples[x, i])
  so <- sample.plots[[s]]
  save_plot(paste0("david-subject-trajectories/frames/", fn, ".png"), so)
}

#' #' # Meta persistent homology
#' #' ## Density level set
#'
#'
#' # density level set -------------------------------------------------------
#'
#'
#' knn.lvls <- sort(unique(V(graf)$mean.kNN))
#' lvl.grafs <- lapply(knn.lvls, function(lvl, g) {
#'   induced_subgraph(g, V(g)$mean.kNN <= lvl)
#' }, g = graf)
#' b0 = sapply(lvl.grafs, function(sg) components(sg)$no)
#' source("run.ends.R")
#' runs <- run.ends(knn.lvls, b0)
#' ggplot(runs, aes(x = start, xend = end, y = value, yend = value)) +
#'   geom_segment() +
#'   labs(x = "density cutoff", y = "b0")
#'
#' xl <- c(min(V(graf)$x) - 2, max(V(graf)$x) + 2)
#' yl <- c(min(V(graf)$y) - 2, max(V(graf)$y) + 2)
#' cl <- sapply(c("min", "max"), function(fn, x) {
#'   do.call(fn, list(x= x))
#'   }, x = V(graf)$mean.kNN)
#' for (i in seq_along(knn.lvls)) {
#'   sg <- lvl.grafs[[i]]
#'   lo <- create_layout(sg, "manual",
#'                       node.positions = data.frame(x = V(sg)$x,
#'                                                   y = V(sg)$y))
#'   ttl <- formatC(knn.lvls[[i]])
#'   p <- plot.mapper(lo, aes_(size = ~size, color = ~mean.kNN),
#'                    list(title = paste("mean kNN =", ttl))) +
#'     coord_cartesian(xlim = xl, ylim = yl) +
#'     scale_color_gradient2(midpoint = mean(V(graf)$mean.kNN),
#'                           high = "blue", mid = "yellow", low = "red",
#'                           limits = cl)
#'   save_plot(paste0("david-level-set/frames/", i, ".png"), p)
#' }
#'
#'
#' # density gradient level set ----------------------------------------------
#'
#' es <- data.table(vnames = attr(E(graf), "vnames"))
#' es[, c("head", "tail") := tstrsplit(vnames, "\\|")]
#' setkey(vertices, name)
#' es[, delta.knn := vertices[tail, mean.knn] - vertices[head, mean.knn]]
#' E(graf)$delta.knn <- es$delta.knn
#' E(graf)$abs.delta <- abs(E(graf)$delta.knn)
#' edge.lvls <- sort(unique(E(graf)$abs.delta))
#' sgs <- lapply(edge.lvls[-1], function(th, graf) {
#'   subgraph.edges(graf, which(E(graf)$abs.delta <= th), delete.vertices = TRUE)
#' }, graf = simplify(graf, remove.multiple = TRUE, edge.attr.comb = "first"))
#' plots <- mapply(function(sg, lvl) {
#'   ttl <- formatC(lvl)
#'   lo <- create_layout(sg, "manual", node.positions = data.frame(
#'     x = V(sg)$x, y = V(sg)$y
#'   ))
#'   plot.mapper(lo, aes_(size = ~size, color = ~mean.kNN),
#'               list(title = paste("thresh =", lvl))) +
#'     coord_cartesian(xlim = xl, ylim = yl) +
#'     scale_color_gradient2(midpoint = mean(V(graf)$mean.kNN),
#'                           high = "blue", mid = "yellow", low = "red",
#'                           limits = cl)
#' }, sg = sgs, lvl = edge.lvls[-1], SIMPLIFY = FALSE)
#' ncomps <- data.table(threshold = edge.lvls[-1],
#'                      b0 = sapply(sgs, count_components))
#' ggplot(ncomps, aes(x = threshold, y = b0)) +
#'   # geom_point() +
#'   geom_line() +
#'   scale_y_log10()
#' for (i in seq_along(plots)) {
#'   save_plot(paste0("david-edge-level-set/frames/", i, ".png"), plots[[i]])
#' }
