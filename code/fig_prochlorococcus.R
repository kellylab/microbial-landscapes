# setup -------------------------------------------------------------------


library(tidyverse)
library(data.table)
library(ggraph)
library(igraph)
library(cowplot)

utils.dir <- "utils/"
for (script in list.files(utils.dir, full.names = TRUE)) source(script)
source("load-prochlorococcus-data.R")
mapper.dir <- "mapper-output/prochlorococcus/"
figs.dir <- "../figures/"
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
    guides(size = FALSE) +
    coord_equal()
}
# set.seed(0)
# lo <- prochlorococcus.mapper$graph %>%
#   activate(nodes) %>%
#   filter(!in.singleton) %>%
#   create_layout("fr", niter = 1000)

# temp ----------------------------------------------------------

plotter(~mean.temp) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill = "C")
temp <- last_plot()
save_plot(paste0(figs.dir, "sup_fig4.pdf"), temp, ncol = 1, nrow = 1,
          base_width = 8)

# depth -------------------------------------------------------------------

plotter(~mean.depth) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_graph(base_family = "Helvetica", base_size = 10) +
  labs(fill = "m") +
  guides(size = FALSE, fill = guide_colorbar(reverse = TRUE))
depth <- last_plot()


# basins ------------------------------------------------------------------

plotter(~as.factor(basin))  +
  labs(fill = "basin")
basins <- last_plot()


# basin series ------------------------------------------------------------

p2basin <- prochlorococcus.v2p[, .N,
                               by = .(site, depth, month, day, cal.month, basin)]
theme_set(theme_cowplot() +
            theme(strip.background = element_rect(size = 0.1),
                  strip.text.x = element_text(size = 10, margin = margin(2, 2, 2, 2)),
                  axis.title = element_text(size = 10),
                  axis.text = element_text(size = 10))
          )
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
      geom_jitter(aes(color = basin), size = 0.1, alpha = 0.1) +
      stat_smooth(aes(fill = basin), data = function(dt) filter(dt, to.smooth)) +
      stat_summary(data = function(dt) filter(dt, !to.smooth),
                   geom = "point", fun.y = mean, size = 0.1) +
      facet_wrap(~ depth, nrow = 4) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_color_hue(drop = FALSE) +
      scale_fill_hue(drop = FALSE) +
      theme(legend.position = "none",
            strip.text.x = element_text(margin = margin(2, 2, 2, 2))) +
      labs(x = "interval (months)", y = "correlation")
  })
pl <- mapply(function(p, nom) p + labs(title = toupper(nom)),
             p = pl, nom = names(pl), SIMPLIFY = FALSE)
plot_grid(plotlist = pl, ncol = 2)
correlation <- last_plot()

# compiled ----------------------------------------------------------------

plot_grid(plot_grid(depth, basins, ncol = 2, labels = c("A", "B")),
          plot_grid(series, correlation, nrow = 2, align = "v", labels = c("C", "D")),
          ncol = 1, rel_heights = c(2, 3))
save_plot(paste0(figs.dir, "fig4.pdf"), last_plot(), ncol = 2,
          base_width = 4, base_height = 10)
