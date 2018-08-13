library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(cowplot)

scripts.dir <- "../r/"
figs.dir <- "../../figures/tda/"
source(paste0(scripts.dir, "load-prochlorococcus-data.R"))
source("k.first.R")
source("vertex.2.points.R")
source("mapper.2.igraph.R")
source("k.first.R")
source("month.2.phase.R")
source("assign.basins.R")
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
                    # mean.phase = mean(phase, na.rm = TRUE),
                    # mean.phase.bats = mean(phase[site == "bats"], na.rm = TRUE),
                    # mean.phase.hot = mean(phase[site == "hot"], na.rm = TRUE),
                    mean.calmonth = mean(cal.month, na.rm = TRUE),
                    median.depth = median(depth, na.rm = TRUE),
                    median.temp = median(temp, na.rm = TRUE),
                    median.sal = median(sal, na.rm = TRUE),
                    mean.knn = mean(knn, na.rm = TRUE)
                    ),
                by = .(vertex, vertex.name)]
graf <- mapper.2.igraph(mpr) %>%
  as_tbl_graph %>%
  left_join(vertices, by = c("name" = "vertex.name"))
graf <- assign.basins(graf, "mean.knn", ignore.singletons = TRUE)
basins <- vertex_attr(graf, "basin") %>%
  as.character %>% unique %>% as.numeric %>% sort(na.last = TRUE) %>%
  as.character
graf <- mutate(graf, basin = factor(basin, levels = basins))
set.seed(1)
lo <- create_layout(graf, "fr", niter = 500)

# paper figure ----------------------------------------------------------------
subplots <- list()
node.maxsize <- 4

#' Composition varies continuously with temperature
subplots$temp <- ggraph(lo) +
  geom_edge_link2(aes(colour = node.mean.temp), show.legend = FALSE) +
  geom_node_point(aes(size = size, color = mean.temp)) +
  scale_color_distiller(palette = "Spectral") +
  scale_edge_color_distiller(palette = "Spectral") +
  scale_size_area(max_size = node.maxsize) +
  labs(color = "C") +
  coord_equal() +
  theme_graph(base_family = "Helvetica", base_size = 10) +
  guides(size = FALSE)

#' Composition varies continuously with depth
subplots$depth <- ggraph(lo) +
  geom_edge_link2(aes(colour = node.mean.depth), show.legend = FALSE) +
  geom_node_point(aes(size = size, color = mean.depth)) +
  scale_color_distiller(palette = "Blues", direction = 1) +
  scale_edge_color_distiller(palette = "Blues", direction = 1) +
  scale_size_area(max_size = node.maxsize) +
  labs(color = "m") +
  coord_equal() +
  theme_graph(base_family = "Helvetica", base_size = 10) +
  guides(size = FALSE)

#' Composition is more stable at low depth and high temperature
subplots$knn <- ggraph(lo) +
  geom_edge_link2(aes(colour = node.mean.knn), show.legend = FALSE) +
  geom_node_point(aes(size = size, color = mean.knn)) +
  scale_color_distiller(palette = "Greys") +
  scale_edge_color_distiller(palette = "Greys") +
  scale_size_area(max_size = node.maxsize) +
  labs(color = "mean\nkNN") +
  coord_equal() +
  theme_graph(base_family = "Helvetica", base_size = 10) +
  guides(size = FALSE)

#' Local kNN minima and basins of attraction
subplots$basins <- ggraph(lo) +
  geom_edge_link0(data = filter(get_edges()(lo), node1.basin != node2.basin |
                                  is.na(node1.basin) | is.na(node2.basin)),
                  colour = "grey") +
  geom_edge_link(aes(colour = node1.basin),
                 data = filter(get_edges()(lo), node1.basin == node2.basin &
                                 !is.na(node1.basin)),
                 show.legend = FALSE) +
  geom_node_point(aes(size = size),
                  data = filter(get_nodes()(lo), is.na(basin)),
                  color = "grey") +
  scale_size_area(max_size = node.maxsize) +
  geom_node_point(aes(size = size, color = basin),
                  data = filter(get_nodes()(lo), !is.na(basin))) +
  geom_node_point(aes(size = size),
                  data = filter(get_nodes()(lo), is.extremum),
                  shape = 21) +
  theme_graph(base_family = "Helvetica") +
  coord_equal() +
  guides(size = FALSE, color = guide_legend(ncol = 2))

#' Dynamics by basin
p2basin <- graf %>% activate(nodes) %>% as.data.table %>%
  merge(v2p, by = "vertex") %>%
  .[, .N, by = .(site, depth, month, day, cal.month, basin)]
theme_set(theme_cowplot(font_size = 10))
pseries <- ggplot(p2basin, aes(x = month, y = basin)) +
  geom_tile(data = function(x) filter(x, is.na(basin)), fill = "grey50") +
  geom_tile(aes(fill = cal.month),
            data = function(x) filter(x, !is.na(basin))) +
  facet_wrap(~ site + depth, nrow = 4) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_y_discrete(limits = c(basins, NA)) +
  scale_fill_gradientn(values = scales::rescale(c(0, 3, 6, 9, 12), c(0, 1)),
                       colours = c("blue", "green", "yellow", "red", "blue")) +
  labs(fill = "calendar\nmonth")
pdistribs <- p2basin[, .(N = sum(N)), by = .(site, depth, cal.month, basin)] %>%
  .[, frac := N / sum(N), by = .(site, depth, cal.month)] %>%
  split(by = "site") %>%
  mapply(function(site, l) {
    ggplot(site, aes(x = cal.month, y = frac)) +
    geom_col(aes(fill = basin), color = "black", size = 0.1, width = 1) +
    facet_wrap(~ depth, nrow = 4) +
    scale_x_continuous(breaks = seq(1, 12, by = 2),
                       labels = as.character(seq(1, 12, by = 2))) +
    theme(legend.position = "none") +
    labs(title = l) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "calendar month", y = "fraction samples")
  }, site = ., l = c("BATS", "HOT"), SIMPLIFY = FALSE) %>%
  plot_grid(plotlist = ., ncol = 2, align = "h",
            axis = "l", vjust = 0.2)
# pad <- plot_grid(pseries, pdistribs, nrow = 2, labels = "AUTO", align = "hv",
#           axis = "lt")

plot_grid(plot_grid(plotlist = subplots, nrow = 2, labels = "AUTO"),
          pdistribs, labels = c("", "E"), nrow = 2,
          rel_heights = c(2, 1))
save_plot(paste0(figs.dir, "paper/fig4.pdf"), last_plot(), base_width = 8,
          nrow = 2)


# other figures -----------------------------------------------------------




#' # Composition is not well-separated by site
site.plots <- mapply(function(d, clr, lo) {
  vs <- unique(d$vertex)
  ggraph(lo) +
    geom_edge_link0(colour = "grey",
                    data = filter(get_edges()(lo),
                                  !(node1.vertex %in% vs) &
                                    !(node2.vertex %in% vs))) +
    geom_edge_link0(data = filter(get_edges()(lo),
                                  node1.vertex %in% vs &
                                    node2.vertex %in% vs),
                    colour = clr) +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size), data = slice(get_nodes()(lo), vs),
                    color = clr)
}, d = split(v2p, by = "site"), clr = c("blue", "red"),
MoreArgs = list(lo = lo), SIMPLIFY = FALSE)
plot_grid(plotlist = site.plots, labels = toupper(names(site.plots)), nrow = 2)
save_plot(paste0(figs.dir, "prochloro-sites.pdf"), last_plot(), base_height = 6,
          nrow = 2)

#' # Composition is not well-separated by salinity
#'
#' BATS is also systematically higher salinity than HOT, so this is skewed.
ggraph(lo) +
  geom_edge_link2(aes(colour = node.mean.sal), show.legend = FALSE) +
  geom_node_point(aes(size = size, color = mean.sal)) +
  scale_color_distiller(palette = "Blues") +
  scale_edge_color_distiller(palette = "Blues") +
  guides(size = FALSE)

#' # Trajectories
sample.xy <- lo %>%
  select(vertex, x, y) %>%
  merge(v2p, by = "vertex") %>%
  merge(as.data.frame(graf), by = "vertex") %>%
  as.data.table %>%
  .[, .(x = mean(x), y = mean(y)), by = .(site, depth, month, year, cal.month)]
setorder(sample.xy, site, depth, month)
sample.jumps <- sample.xy %>%
  split(by = c("site", "depth")) %>%
  lapply(function(d, lo) {
    setorder(d, month)
    ggraph(lo) +
      geom_edge_link0(colour = "grey") +
      geom_node_point(aes(size = size), colour = "grey") +
      geom_path(aes(x = x, y = y, color = cal.month / 12), data = d,
                arrow = arrow(angle = 15, type = "open",
                              length = unit(10, "points"))
                ) +
      geom_point(aes(x = x, y = y, color = cal.month / 12), data = d) +
      season_gradient +
      guides(size = FALSE) +
      labs(color = "month i")
  }, lo = lo)
sample.jumps <- unique(prochlorococcus[, .(site, depth)]) %>%
  setorder(site, depth) %>%
  .[, paste(site, depth, sep = ".")] %>%
  sample.jumps[.]
save_plot(paste0(figs.dir, "prochloroccocus-steps.pdf"),
  plot_grid(plotlist = sample.jumps, labels = names(sample.jumps),
            nrow = 4, ncol = 6),
  nrow = 4, ncol = 6, base_height = 4
  )
