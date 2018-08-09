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
source("assign.basins.R")
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
k <- round(nrow(samples) / 10)
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

# mapper call -------------------------------------------------------------

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

#' Find basins of attraction:
graf <- assign.basins(graf, "mean.knn", ignore.singletons = TRUE)
graf <- graf %>% 
  activate(nodes) %>% 
  mutate(basin = factor(basin,
                        levels = as.character(sort(as.numeric(unique(basin))),
                                              na.last = TRUE)))

#' ## Plots

# paper figure ----------------------------------------------------------------


#' Fraction of samples in each vertex belonging to each subject
set.seed(0)
lo <- create_layout(graf, "fr")
subplots <- list()
subplots$fsubject <- ggraph(lo) +
  geom_edge_link2(aes(colour = node.subject), show.legend = FALSE) +
  geom_node_point(aes(size = size, color = subject)) +
  labs(color = "fraction A") +
  scale_color_distiller(palette = "Spectral") +
  scale_edge_colour_distiller(palette = "Spectral") +
  theme_graph(base_family = "Helvetica") +
  coord_equal() +
  guides(size = FALSE)
# save_plot(paste0(figs.dir, "david-fsubject.pdf"), last_plot(), base_height = 6)

#' kNN
subplots$knn <- ggraph(lo) +
  geom_edge_link2(aes(colour = node.mean.knn), show.legend = FALSE) +
  geom_node_point(aes(size = size, color = mean.knn)) +
  labs(color = "mean\nkNN") +
  scale_color_distiller(palette = "Blues") +
  scale_edge_colour_distiller(palette = "Blues") +
  theme_graph(base_family = "Helvetica") +
  coord_equal() +
  guides(size = FALSE)

#' ### Basins of attraction
subplots$basins <- ggraph(lo) +
  geom_edge_link0(data = filter(get_edges()(lo), node1.basin != node2.basin |
                                  is.na(node1.basin) | is.na(node2.basin)),
                  colour = "grey50") +
  geom_edge_link0(aes(colour = node1.basin),
                  data = filter(get_edges()(lo), node1.basin == node2.basin),
                  show.legend = FALSE) +
  geom_node_point(aes(size = size),
                  data = filter(get_nodes()(lo), is.na(basin)), color = "grey50") +
  geom_node_point(aes(size = size, color = basin),
                  data = filter(get_nodes()(lo), !is.na(basin))) +
  geom_node_point(aes(size = size),
                  data = filter(get_nodes()(lo), is.extremum),
                  shape = 21) +
  coord_equal() +
  theme_graph(base_family = "Helvetica") +
  theme(legend.position = "none") +
  guides(size = FALSE)
# save_plot(paste0(figs.dir, "david-basins.pdf"), last_plot(), base_height = 6)

#' ## Subject trajectories by basin
sample.basins <- graf %>%
  activate(nodes) %>%
  select(-subject) %>%
  as.data.table %>%
  merge(v2p, by = "vertex") %>%
  .[, .N, by = .(subject, day, point, event, basin)]  
sample.basins[, i := rank(basin, ties.method = "first"), by = point]
theme_set(theme_cowplot(font_size = 8))
basins <- sample.basins$basin %>% 
  unique %>% 
  sort(na.last = TRUE)
sample.basins[, i := rank(basin, ties.method = "first"), by = point]
pseries <- sample.basins %>% 
  ggplot(aes(x = day, y = basin)) +
  geom_point(data = function(dt) filter(dt, is.na(basin)),
             color = "grey50", shape = 21, size = 1) +
  geom_point(aes(color = event),
             data = function(dt) filter(dt, !is.na(basin)), size = 1) +
  # geom_tile(aes(fill = basin)) +
  scale_y_discrete(limits = basins) +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom") +
  guides(color = guide_legend(direction = "horizontal", nrow = 2)) +
  background_grid(major = "y") +
  facet_grid(subject ~ .)

#' Distribution across basins during different events
pdistribs <- sample.basins %>% 
  group_by(subject, event) %>%
  mutate(frac = N / sum(N)) %>%
  group_by(subject, event, basin) %>% 
  summarize(frac = sum(frac)) %>% 
  ggplot(aes(x = event)) +
  geom_col(aes(y = frac, group = basin, fill = basin),
           position = "stack") +
  coord_flip() +
  labs(y = "fraction samples") +
  guides(fill = guide_legend(direction = "horizontal", nrow = 1,
                             label.position = "bottom")) +
  facet_grid(subject ~ ., scales = "free")
basin.legend <- get_legend(pdistribs)
pdistribs <- pdistribs + theme(legend.position = "none")
subplots$basin.distribs <- plot_grid(pseries, pdistribs, 
                                     nrow = 2, #align = "hv", axis = "lt",
                                     labels = c("D", "E"), 
                                     rel_heights = c(1.5, 1))
plot_grid(plot_grid(plotlist = subplots, labels = c("A", "B", "C", NA), nrow = 2), 
          basin.legend, nrow = 2, rel_heights = c(5, 1))
save_plot(paste0(figs.dir, "paper/fig3.pdf"), last_plot(), nrow = 2, ncol = 2)


# other figures -----------------------------------------------------------



#' JSDs of distribution across vertices during different events
v2p[, frac := 1 / .N, by = point]
event.vdist <- v2p[, .(N = sum(frac)), by = .(subject, event, vertex)] #%>%
event.densities <- event.vdist %>%
  split(by = c("subject", "event")) %>%
  lapply(function(event, lo) {
    nodes <- left_join(get_nodes()(lo), event, by = "vertex")
    ggraph(lo) +
      geom_edge_link0(colour = "grey") +
      geom_node_point(data = filter(nodes, is.na(N)), color = "grey") +
      geom_node_point(aes(color = N), data = filter(nodes, !is.na(N))) +
      scale_color_distiller(palette = "Spectral") +
      scale_alpha(na.value = 0) +
      coord_equal() +
      guides(size = FALSE)

  }, lo = lo)
plot_grid(plotlist = event.densities, labels = names(event.densities),
          ncol = 5, align = "hv")
save_plot(paste0(figs.dir, "david-event-densities.pdf"),
          last_plot(), ncol = 5, base_height = 6,
          base_aspect_ratio = 0.7)
event.vdist <- event.vdist %>%
  dcast(subject + event ~ vertex, value.var = "N", fill = 0)
se <- paste(event.vdist$subject, event.vdist$event, sep = ".")
event.vdist <- as.matrix(event.vdist[, -c(1, 2)])
rownames(event.vdist) <- se
vdist.jsd <- JSD(event.vdist, est.prob = "empirical")
rownames(vdist.jsd) <- se
colnames(vdist.jsd) <- se
melt(vdist.jsd, varnames = c("event.i", "event.j"), value.name = "jsd") %>%
  ggplot(aes(x = event.i, y = event.j)) +
  geom_tile(aes(fill = jsd)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_distiller(palette = "Blues", direction = -1) +
  coord_equal()

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
