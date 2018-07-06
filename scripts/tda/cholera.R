library(data.table)
library(TDAmapper)
library(philentropy)
library(cowplot)
library(tidyr)
library(igraph)
library(ggraph)
library(tidygraph)

util.dir <- "../r/"

source(paste0(util.dir, "load_cholera_data.R"))
source("vertex.2.points.R")
source("dist2knn.R")
source("mapper.2.igraph.R")
source("plot.mapper.R")
source("sample.subgraphs.R")

# load data ---------------------------------------------------------------

gordon[, freq := count / sum(count), by = sample]
gordon.samples <- unique(gordon[, .(sample, subject, diagnosis, id, hour)])
gordon.samples[, idx := frank(hour), by = subject]
gordon.samples[, progression := idx / max(idx), by = subject]
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

#' # K-nearest neighbor density

# knn ---------------------------------------------------------------------


k <- 10
knn <- dist2knn(js.dist, k = k)
gordon.samples[, knn := knn[sample]]
ggplot(gordon.samples, aes(x = progression, y = knn)) +
  geom_smooth() +
  geom_point(aes(color = diagnosis)) +
  facet_wrap(~ subject)

#' Diarrhea samples only showing disease duration:
setkey(gordon.samples, diagnosis)
ggplot(gordon.samples["diarrhea"], aes(x = hour, y = knn)) +
  geom_smooth() +
  geom_point() +
  facet_wrap(~ subject)

# mds --------------------------------------------------------

mds2 <- cmdscale(js.dist, 2, eig = TRUE)
mds2$GOF # lossy

#' MDS sketch suggests existence of at least 2 clusters.
plot(mds2$points)

#' MDS rank evens out density but preserves some of the patchiness
rk.mds <- apply(mds2$points, 2, rank)
plot(rk.mds)

po <- 70
ni <- c(10, 10)
mpr <- mapper2D(js.dist, filter_values = list(rk.mds[,1], rk.mds[,2]),
                percent_overlap = po,
                num_intervals = ni)
v2p <- vertex.2.points(mpr$points_in_vertex)
v2p <- merge(v2p, gordon.samples, by.x = "point.name", by.y = "sample")
vertices <- v2p[, .(mean.knn = mean(knn),
                    f.state = sum(diagnosis == "diarrhea") / .N,
                    mean.t = mean(hour),
                    size = .N),
                by = .(vertex, vertex.name)]
setorder(vertices, vertex)
graf <- mapper.2.igraph(mpr)
graf <- graf %>%
  as_tbl_graph %>%
  activate(nodes) %>%
  left_join(vertices, by = c("name" = "vertex.name"))

# draw graphs -------------------------------------------------------------

set.seed(1)
lo <- create_layout(graf, "fr", niter = 1000)
# color by fraction samples marked diarrhea
plot.mapper(lo, aes_(size = ~size, color = ~f.state),
            list(size = "# samples", color = "f diarrhea")) +
  scale_color_distiller(palette = "Spectral")
# color vertices by mean knn density
plot.mapper(lo, aes_(size = ~size, color = ~mean.knn),
            list(size = "# samples", color = "mean kNN")) +
  scale_color_distiller(palette = "Blues")


#' # Subject trajectories
# subject trajectories ----------------------------------------------------
graf <- graf %>%
  activate(nodes) %>%
  mutate(x = lo$x, y = lo$y)
subject.phases <- v2p %>%
  split(by = "subject") %>%
  lapply(function(dt, graf) {
    max.i <- max(dt$idx)
    v <- dt[, .(progression = mean(progression)), by = .(vertex, vertex.name)]
    sg <- graf %>%
      slice(v$vertex)
    sg %>%
      activate(nodes) %>%
      left_join(v, by = c("name" = "vertex.name"))

  }, graf = graf)
subject.phase.plots <- mapply(function(g, subj, global) {
  d <- as.data.frame(g)
  plot.mapper(global, aes_(size = ~size), labs(title = subj), color = "grey") +
    geom_point(aes(x = x, y = y, size = size, color = progression), data = d) +
    scale_color_distiller(palette = "Spectral", direction = 1)
}, g = subject.phases, subj = names(subject.phases),
MoreArgs = list(global = lo), SIMPLIFY = FALSE)
for (i in seq_along(subject.phase.plots)) {
  fn <- paste0("cholera-subject-trajectories/phases/",
               names(subject.phase.plots)[i], ".png")
  save_plot(fn, subject.phase.plots[[i]], base_height = 6)
}

#' ## "Top" arc progressions
cl1 <- c("A", "F", "G")
plot_grid(plotlist = subject.phase.plots[cl1], labels = cl1)

#' ## Other
cl.other <- c("B", "C", "D", "E")
plot_grid(plotlist = subject.phase.plots[cl.other], labels = cl.other)

# frames ------------------------------------------------------------------


sample.sgs <- v2p %>%
  split(by = "point.name") %>%
  lapply(function(s, graf) {
    slice(graf, s$vertex)
  }, graf = graf)
sample.frames <- mapply(function(sg, sample, state, graf) {
  if (state == "diarrhea") {
    color <- "red"
  } else {
    color <- "blue"
  }
  d <- as.data.frame(sg)
  plot.mapper(graf, aes_(size = ~size), labs(title = sample), color = "grey") +
    geom_point(aes(x = x, y = y, size = size), data = d, color = color)
}, sg = sample.sgs[gordon.samples$sample], sample = gordon.samples$sample,
state = gordon.samples$diagnosis,
MoreArgs = list(graf = lo), SIMPLIFY = FALSE)
setorder(gordon.samples, subject, idx)
for (i in seq(nrow(gordon.samples))) {
  fn <- paste0("cholera-subject-trajectories/frames/",
               gordon.samples[i, subject], gordon.samples[i, idx],
               ".png")
  save_plot(fn, sample.frames[[gordon.samples[i, sample]]])
}

#' #' # 'Meta' persistent homology on the Mapper graph
#' #' ## b0 vs density threshold
# # b0 vs density threshold ----------------------------------------------
#
# ls <- sort(unique(V(g1)$mean.knn))
# ls.gs <- lapply(ls, function(x, g) induced_subgraph(g, V(g)$mean.knn <= x),
#                 g = g1)
# comps <- lapply(ls.gs, components)
# # ncomps <- sapply(comps, function(g) sum(g$csize > 1)) # filter singletons
# ncomps <- sapply(comps, function(g) g$no)
# ggplot(data.frame(threshold = ls, b0 = ncomps), aes(x = threshold, y = b0)) +
#   geom_line() + geom_point()
#
#' #' Get the longest interval in density threshold where b0 is constant.
#' #' First get the number of intervals in the threshold for each run where b0 is
#' #' constant:
# b0.runs <- rle(ncomps)
# b0.runs <- data.table(length = b0.runs$lengths, b0 = b0.runs$values)
# #' Cumulatively sum over intervals in the threshold value to get the end
# #' threshold value of each run:
# b0.runs[, end := sapply(cumsum(length), function(n, th) {
#   th[n]
# }, th = ls)]
# #' Shift the end values down and prepend 0 for the start threshold value of
# #' each run:
# b0.runs[, start := c(0, end[-.N])]
# #' Get persistence length of each run in the dimension of the threshold.
# b0.runs[, persistence := end - start]
# #' Get rid of the run between 0 and the first point.
# b0.runs <- b0.runs[-1,]
# #' Find the threshold at the end of the run with maximum persistence
# thc <- b0.runs[b0 > 1, end[persistence == max(persistence)]]
#
# #' Get the vertices in that level set:
# rep.vs <- which(V(g1)$mean.knn <= thc)
# rep.xy <- lo[rep.vs, c("x", "y")]
# rep.g <- induced_subgraph(g1, rep.vs)
# V(rep.g)$state <- as.character(comps[[which(ls == thc)]]$membership)
# rep.lo <- create_layout(rep.g, "manual", node.positions = rep.xy)
# plot.mapper(rep.lo, aes_(size = ~size, color = ~state), list(size = "samples"))
#
# #' ## Video of graph expansion as vertices are added from least to most dense
# # global density level set video frames -----------------------------------
# lvlset.grafs <- lapply(ls, function(th, graf, layout) {
#   vids <- which(V(graf)$mean.knn <= th)
#   sg <- induced_subgraph(graf, vids)
#   xl <- c(min(V(graf)$x) - 1, max(V(graf)$x) + 1)
#   yl <- c(min(V(graf)$y) - 1, max(V(graf)$y) + 1)
#   cl <- sapply(c("min", "max"), function(fn, x) {
#     do.call(fn, list(x= x))
#     }, x = V(graf)$mean.knn)
#   ggraph(sg, "manual", node.positions = data.frame(x = V(sg)$x, y = V(sg)$y)) +
#     geom_edge_link0() +
#     geom_node_point(aes(color = mean.knn, size = size)) +
#     coord_cartesian(xlim = xl, ylim = yl) +
#     scale_color_gradient2(midpoint = median(V(graf)$mean.knn),
#                           high = "blue", mid = "yellow", low = "red",
#                           limits = cl) +
#     theme_graph()
# }, graf = g1, layout = lo)
# for (i in seq_along(ls)) {
#   fn <- paste0("level-set/frames/", i, ".png")
#   save_plot(fn, lvlset.grafs[[i]])
# }
#
