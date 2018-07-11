library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(cowplot)
library(combinat)

scripts.dir <- "../r/"
figs.dir <- "../../figures/tda/"
source(paste0(scripts.dir, "load-prochlorococcus-data.R"))
source("k.first.R")
source("vertex.2.points.R")
source("mapper.2.igraph.R")
source("plot.mapper.R")
source("k.first.R")
source("month.2.phase.R")
source("layout.tbl.graph.R")
source("as.layout.manual.R")
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
set.seed(1)
graf <- layout.tbl.graph(graf, "fr", niter = 500)

#' # Composition varies continuously with temperature
plot.mapper(as.layout.manual(graf), aes_(size = ~size, color = ~mean.temp)) +
  scale_color_distiller(palette = "Spectral")

#' # Composition varies continuously with depth
plot.mapper(as.layout.manual(graf), aes_(size = ~size, color = ~mean.depth)) +
  scale_color_distiller(palette = "Blues", direction = 1)
save_plot(paste0(figs.dir, "prochloro-depth.pdf"), last_plot(), base_height = 6)

#' # Composition is more stable at low depth and high temperature
plot.mapper(as.layout.manual(graf), aes_(size = ~size, color = ~mean.knn)) +
  scale_color_distiller(palette = "Spectral")

#' # Composition is not well-separated by site
site.plots <- mapply(function(d, clr, graf) {
  lo <- as.layout.manual(graf)
  sg <- slice(lo, unique(d$vertex))
  df <- as.data.frame(sg)
  plot.mapper(lo, aes_(size = ~size), NULL, color = "grey") +
    geom_point(aes(x = x, y = y, size = size), data = df, color = clr)
}, d = split(v2p, by = "site"), clr = c("blue", "red"),
MoreArgs = list(graf = graf), SIMPLIFY = FALSE)
plot_grid(plotlist = site.plots, labels = toupper(names(site.plots)), nrow = 2)
save_plot(paste0(figs.dir, "prochloro-sites.pdf"), last_plot(), base_height = 6,
          nrow = 2)

#' # Composition is not well-separated by salinity
#'
#' BATS is also systematically higher salinity than HOT, so this is skewed.
plot.mapper(as.layout.manual(graf), aes_(size = ~size, color = ~mean.sal)) +
  scale_color_distiller(palette = "Spectral")

#' # Monthwise gradients per site + depth
season_gradient <- scale_color_gradientn(
  values = scales::rescale(c(1, 3, 6, 9, 12), c(0, 1)),
  colors = c("blue", "green", "yellow", "red", "blue"),
  labels = c("1", "3", "6", "9", "12"),
  limits = c(0, 1))
depth.phase <- v2p %>%
  split(by = c("site", "depth"), flatten = FALSE) %>%
  lapply(function(d) {
    lapply(d, function(d, graf) {
      d2 <- d[, .(mean.phase = mean(cal.month) /  12), by = vertex]
      sg <- as.data.frame(slice(graf, d2$vertex))
      sg <- merge(sg, d2, by = "vertex")
      # browser()
      plot.mapper(as.layout.manual(graf), aes_(size = ~size), NULL,
                  color = "grey") +
       geom_point(aes(x = x, y = y, size = size, color = mean.phase), data = sg) +
        guides(size = FALSE) +
        season_gradient +
        labs(color = "mean month")
      }, graf = graf)
  })
for (si in seq_along(depth.phase)) {
  depths <- names(depth.phase[[si]]) %>%
    as.numeric %>%
    sort %>%
    as.character
  depth.phase[[si]] <- depth.phase[[si]][depths]
  pp <- plot_grid(plotlist = depth.phase[[si]], labels = depths)
  fn <- paste0(figs.dir, "prochloro-phase-", names(depth.phase)[si], ".pdf")
  save_plot(fn, pp, nrow = 4, ncol = 4, base_aspect_ratio = 1.5)
}

ncomponents <- v2p %>%
  split(by = c("site", "depth", "cal.month")) %>%
  lapply(function(d, graf) slice(graf, unique(d$vertex)), graf = graf) %>%
  lapply(count_components) %>%
  as.data.frame %>%
  melt(id.vars = NULL, variable.name = "bin", value.name = "ncomps") %>%
  as.data.table
ncomponents[, c("site", "depth", "cal.month") := tstrsplit(bin, "\\.", type.convert = TRUE)]
ggplot(ncomponents, aes(x = cal.month, y = ncomps)) +
  geom_hline(yintercept = 1) +
  geom_point(aes(color = depth)) +
  scale_color_gradient(low = "cyan2", high = "darkblue") +
  facet_wrap(~ site + depth)
