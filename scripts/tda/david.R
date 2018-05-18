library(philentropy)
library(TDAmapper)
library(data.table)
library(tidyverse)
library(ggraph)
library(igraph)
library(cowplot)

# setup -------------------------------------------------------------------


scripts.dir <- "../r/"
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
source("dist2knn.R")
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
#' 2D MDS shows much more variance in dimension 1 than dimension 2, resulting in
#' large empty spaces.


# mds ------------------------------------------------------------------


mds2d <- as.data.table(cmdscale(dist.mat, 2), keep.rownames = "sample")
setnames(mds2d, c("V1", "V2"), c("mds1", "mds2"))
mds2d[, subject := tstrsplit(sample, "_")[[1]]]
ggplot(mds2d, aes(x = mds1, y = mds2)) +
  geom_point(aes(color = subject))

#' With log kNN density we can spread out a bit.

mds1d <- cmdscale(dist.mat, 1)
samples[, mds1 := mds1d[sample,1]]
ggplot(samples, aes(x = mds1, y = kNN)) +
  geom_point(aes(color = subject)) +
  scale_y_log10()

#' Or with the time:
samples[, scaled.time := day / max(day), by = subject]
ggplot(samples, aes(x = mds1, y = scaled.time)) +
  geom_point(aes(color = subject))

#' # Mapper
#' ## Filter by 2D MDS

# mds filter -------------------------------------------------------------


mpr <- mapper2D(dist.mat, list(mds2d$mds1, mds2d$mds2),
                num_intervals = c(100, 10),
                percent_overlap = 80)
mpr$points_in_vertex <- lapply(mpr$points_in_vertex, function(v, rn) {
  names(v) <- rn[v]
  v
}, rn = rownames(dist.mat))


# mapper output ----------------------------------------------------


# source("mapper.adj.R")
# mpr$adjacency <- mapper.adj(mpr$points_in_vertex) # correct adj matrix

#' Format Mapper output:
source("vertex.2.points.R")
v2p <- vertex.2.points(mpr$points_in_vertex) # map vertices to samples
v2p <- merge(v2p, samples, by.x = "point.name", by.y = "sample")
source("mapper.2.igraph.R")
graf <- mapper.2.igraph(mpr)
V(graf)$subject <- sapply(mpr$points_in_vertex, function(pts, dt) {
  pnames <- names(pts)
  setkey(dt, sample)
  dt[pnames, sum(subject == "A") / .N]
}, dt = samples)


#' ## Plots

source("plot.mapper.R")
#' Fraction of samples in each vertex belonging to each subject:
set.seed(0)
lo <- create_layout(graf, "fr")
# store xy positions in graph for later
V(graf)$x <- lo$x
V(graf)$y <- lo$y
plot.mapper(lo, aes_(size = ~size, color = ~subject),
            list(color = "fraction A")) +
  scale_color_gradient2(midpoint = 0.5, mid = "yellow")

#' ### Density
setkey(samples, sample)
V(graf)$mean.kNN <- sapply(V(graf)$name, function(v, dt) {
  setkey(dt, vertex.name)
  dt[v, mean(kNN)]
}, dt = v2p)
plot.mapper(create_layout(graf, "manual",
                          node.positions = data.frame(x = V(graf)$x,
                                                      y = V(graf)$y)),
            aes_(size = ~size, color = ~mean.kNN)) +
  scale_color_distiller(palette = "Spectral", trans = "log10")

#' ### Events
source(paste0(scripts.dir, "Mode.R"))
subj.grafs <- lapply(c("A", "B"), function(subj, graf, v2p) {
  setkey(v2p, subject)
  g <- induced_subgraph(graf, v2p[subj, unique(vertex)])
  setkey(v2p, vertex.name)
  V(g)$event <- sapply(V(g)$name, function(v, subj, v2p) {
    setkey(v2p, subject, vertex.name)
    paste(Mode(v2p[.(subj, v), event]), collapse = "/")
  }, subj = subj, v2p = v2p)
  g
}, graf = graf, v2p = v2p)
subj.plots <- mapply(function(graf, subj) {
  plot.mapper(create_layout(graf, "manual",
                            node.positions = data.frame(x = V(graf)$x,
                                                        y = V(graf)$y)),
              aes_(size = ~size, color = ~event),
              list(title = subj)) +
    theme(aspect.ratio = 1)
}, graf = subj.grafs, subj = c("A", "B"), SIMPLIFY = FALSE)
plot_grid(plotlist = subj.plots)


#' ### Eric's trip
#'
# eric’s frames -----------------------------------------------------------


source("sample.subgraphs.R")
sample.sgrafs <- sample.subgraphs(
  v2p[, .(samples = point.name, vertices = vertex)],
  graf)
# setkey(v2p, subject, day)
# setkey(events, subject, day)
# v2p <- events[v2p]
  setkey(samples, subject)
b.events <- samples["B", .(sample, event)]
sample.overlays <- mapply(function(s, sg, e, lo) {
  if (!is.null(sg)) {
    if (grepl("pre", e)) {
      clr <- "blue"
    } else if (grepl("post", e)) {
      clr <- "yellow"
    } else {
      clr <- "red"
    }
    ttl <- paste(s, e, sep = "; ")
    g <- plot.mapper(lo, aes_(size = ~size), list(title = ttl),
                     color = "grey50")
    g + geom_point(aes(x = x, y = y, size = size),
                   data = as.data.frame(vertex_attr(sg)),
                   color = clr)



  } else {
    NULL
  }
}, s = b.events$sample, sg = sample.sgrafs[b.events$sample],
e = b.events$event,
MoreArgs = list(lo = lo), SIMPLIFY = FALSE)
setorder(samples, day)
setkey(samples, subject)
sample.overlays <- sample.overlays[samples["B", sample]]
sample.overlays <- sample.overlays[sapply(sample.overlays,
                                          function(x) !is.null(x))]

# write -------------------------------------------------------------------


for (i in seq_along(sample.overlays)) {
  so <- sample.overlays[[i]]
  save_plot(paste0("david-subject-trajectories/frames/B", i, ".png"), so)
}

#' # Meta persistent homology
#' ## Density level set


# density level set -------------------------------------------------------


knn.lvls <- sort(unique(V(graf)$mean.kNN))
lvl.grafs <- lapply(knn.lvls, function(lvl, g) {
  induced_subgraph(g, V(g)$mean.kNN <= lvl)
}, g = graf)
b0 = sapply(lvl.grafs, function(sg) components(sg)$no)
source("run.ends.R")
runs <- run.ends(knn.lvls, b0)
ggplot(runs, aes(x = start, xend = end, y = value, yend = value)) +
  geom_segment() +
  labs(x = "density cutoff", y = "b0")

xl <- c(min(V(graf)$x) - 2, max(V(graf)$x) + 2)
yl <- c(min(V(graf)$y) - 2, max(V(graf)$y) + 2)
cl <- sapply(c("min", "max"), function(fn, x) {
  do.call(fn, list(x= x))
  }, x = V(graf)$mean.kNN)
for (i in seq_along(knn.lvls)) {
  sg <- lvl.grafs[[i]]
  lo <- create_layout(sg, "manual",
                      node.positions = data.frame(x = V(sg)$x,
                                                  y = V(sg)$y))
  ttl <- formatC(knn.lvls[[i]])
  p <- plot.mapper(lo, aes_(size = ~size, color = ~mean.kNN),
                   list(title = paste("mean kNN =", ttl))) +
    coord_cartesian(xlim = xl, ylim = yl) +
    scale_color_gradient2(midpoint = mean(V(graf)$mean.kNN),
                          high = "blue", mid = "yellow", low = "red",
                          limits = cl)
  save_plot(paste0("david-level-set/frames/", i, ".png"), p)
}
