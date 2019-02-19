# setup -------------------------------------------------------------------


library(tidyverse)
library(data.table)
library(ggraph)
library(igraph)
library(cowplot)

utils.dir <- "utils/"
for (script in list.files(utils.dir, full.names = TRUE)) source(script)
loader.script <- "../r/load-prochlorococcus-data.R"
source(loader.script)
mapper.dir <- "../../output/mapper/prochlorococcus/"
prochlorococcus.mapper <- read.mapper.graph(mapper.dir)
prochlorococcus.mapper$graph <- prochlorococcus.mapper$graph %>%
  activate(nodes) %>%
  mutate(basin = as.factor(basin))
prochlorococcus.v2p <- merge.mpr.samples(prochlorococcus.mapper,
                                         prochlorococcus.samples)
plotter <- function(v) {
  plot.mapper.graph(prochlorococcus.mapper$graph,
                    node = geom_node_point(aes_(size = ~size, fill = v),
                                           shape = 21),
                    seed = 0,
                    exclude.singletons = TRUE) +
    guides(size = FALSE)
}
# set.seed(0)
# lo <- prochlorococcus.mapper$graph %>%
#   activate(nodes) %>%
#   filter(!in.singleton) %>%
#   create_layout("fr", niter = 1000)

# temp ----------------------------------------------------------


# theme_set(theme_graph(base_family = "Helvetica", base_size = 10))
#' Composition varies continuously with temperature
plotter(~mean.temp) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill = "C")

ggraph(lo) +
  geom_edge_link0() +
  geom_node_point(aes(size = size, fill = mean.temp), shape = 21) +
  scale_fill_distiller(palette = "Spectral") +
  theme_graph(base_family = "Helvetica", base_size = 10) +
  labs(fill = "C") +
  coord_equal() +
  guides(size = FALSE)
temp <- last_plot()

# depth -------------------------------------------------------------------



#' Composition varies continuously with depth
ggraph(lo) +
  geom_edge_link0() +
  geom_node_point(aes(size = size, fill = mean.depth), shape = 21) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_graph(base_family = "Helvetica", base_size = 10) +
  labs(fill = "m") +
  coord_equal() +
  guides(size = FALSE, fill = guide_colorbar(reverse = TRUE))
depth <- last_plot()


# basins ------------------------------------------------------------------


ggraph(lo) +
  geom_edge_link0() +
  geom_node_point(aes(size = size, fill = basin), shape = 21) +
  theme_graph(base_family = "Helvetica", base_size = 10) +
  coord_equal() +
  guides(size = FALSE, fill = guide_legend(ncol = 2))
basins <- last_plot()


# basin series ------------------------------------------------------------

p2basin <- prochlorococcus.v2p[, .N,
                               by = .(site, depth, month, day, cal.month, basin)]
theme_set(theme_cowplot())
plot.series <- function(dt) {
  jans <- dt[cal.month == 1, unique(month)]
  ggplot(dt, aes(x = month, y = basin)) +
    geom_tile(aes(fill = basin)) +
    geom_vline(xintercept = jans, linetype = "dotted") +
    facet_wrap(~ depth, nrow = 4) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "none")
}
pl <- p2basin %>%
  split(by = "site") %>%
  lapply(plot.series)
pl <- mapply(function(p, n) p + labs(title = toupper(n)), p = pl, n = names(pl),
             SIMPLIFY = FALSE)
plot_grid(plotlist = pl)
series <- last_plot()

# distribs ----------------------------------------------------------------

p2basin %>%
  .[, .(N = sum(N)), by = .(site, depth, basin, cal.month)] %>%
  .[, frac := N / sum(N), by = .(site, depth, cal.month)] %>%
  split(by = "site") %>%
  mapply(function(site, l) {
    ggplot(site, aes(x = cal.month, y = frac)) +
      geom_col(aes(fill = basin), color = "black", size = 0.1, width = 1) +
      facet_wrap(~ depth, nrow = 4) +
      scale_x_continuous(breaks = seq(1, 12, by = 2),
                         labels = as.character(seq(1, 12, by = 2))) +
      labs(title = l) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      theme(legend.position = "none") +
      labs(x = "calendar month", y = "fraction samples")
  }, site = ., l = c("BATS", "HOT"), SIMPLIFY = FALSE) %>%
  plot_grid(plotlist = ., ncol = 2, align = "h",
            axis = "l", vjust = 0.2)
distribs <- last_plot()


# correlation -------------------------------------------------------------

prochlorococcus.persistence <- prochlorococcus.v2p %>%
  split(by = c("site", "depth")) %>%
  lapply(function(df) basin.persistence(df$basin, df$month, scale = TRUE)) %>%
  rbindlist(idcol = "site.depth") %>%
  .[, c("site", "depth") := tstrsplit(site.depth, "\\.")] %>%
  .[, site.depth := NULL] %>%
  .[, basin := as.factor(as.numeric(basin))] %>%
  .[, depth := as.numeric(depth)]
# filter out depth-basin combos that have < 10 unique delta.t values
# these mess up the smoothing calculation for the entire facet!
prochlorococcus.persistence[, to.smooth := uniqueN(delta.t) >= 10,
                            by = .(depth, basin)]
pl <- prochlorococcus.persistence[!is.na(basin)] %>%
  split(by = "site") %>%
  lapply(function(df) {
    ggplot(df, aes(x = delta.t, y = f, color = basin)) +
      stat_smooth(aes(fill = basin), data = function(dt) filter(dt, to.smooth)) +
      stat_summary(data = function(dt) filter(dt, !to.smooth),
                   geom = "point", fun.y = mean, size = 0.1) +
      facet_wrap(~ depth) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_color_hue(drop = FALSE) +
      scale_fill_hue(drop = FALSE) +
      labs(x = "interval (months)", y = "correlation")
  })
pl <- mapply(function(p, nom) p + labs(title = toupper(nom)),
             p = pl, nom = names(pl), SIMPLIFY = FALSE)
plot_grid(plot_grid(pl[[1]] + theme(legend.position = "none"),
                    pl[[2]] + theme(legend.position = "none"),
                    ncol = 2))
correlation <- last_plot()

# compiled ----------------------------------------------------------------

plot_grid(plot_grid(temp, depth, nrow = 1, labels = "AUTO"),
          basins, labels = c("", "C"),
          ncol = 1)
save_plot("../../figures/tda/paper/fig6.pdf", last_plot(),
          ncol = 1,
          base_width = 8, base_height = 8)
plot_grid(series, distribs, correlation,
          ncol = 1, labels = "AUTO")
save_plot("../../figures/tda/paper/fig7.pdf", last_plot(), ncol = 1,
          base_width = 8, base_height = 12)
