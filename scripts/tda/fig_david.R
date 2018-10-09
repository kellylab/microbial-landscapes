# setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(cowplot)

loader.script <- "../r/load_david_data.R"
mapper.dir <- "../../output/mapper/david/"
utils.dir <- "utils/"
for (script in list.files(utils.dir, full.names = TRUE)) source(script)
source(loader.script)
david.mapper <- read.mapper.graph(mapper.dir)
david.v2p <- merge.mpr.samples(david.mapper, david.samples)

# frac subject ------------------------------------------------------------


set.seed(13)
lo <- david.mapper$graph %>% # mapper graph with no outlier
  activate(nodes) %>%
  filter(!in.singleton) %>%
  create_layout("fr", niter = 1000)
ggraph(lo) +
  geom_edge_link0() +
  geom_node_point(aes(fill = f.subject, size = size),
                  data = function(df) filter(df, !in.singleton), shape = 21) +
  labs(fill = "subject") +
  scale_fill_distiller(palette = "Spectral", breaks = c(0, 0.5, 1),
                       labels = c("all B", "", "all A")) +
  scale_size_area(max_size = 4) +
  guides(size = FALSE) +
  theme_graph(base_family = "Helvetica") +
  coord_equal()
fsubject <- last_plot()


# events ------------------------------------------------------------------


setkey(david.v2p, subject, event)
title.size <- 10
base.size <- 8
theme_set(theme_graph(base_family = "Helvetica", title_size = title.size,
                      base_size = base.size) +  theme(legend.position = "none"))
plot_grid(
  ggraph(lo) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size),
                    data = function(df) {
                      filter(df, vertex %in% david.v2p[.("A", "US (pre)")]$vertex)
                    },
                    fill = "blue", shape = 21) +
    labs(title = "A, pre-travel") +
    scale_size_area(max_size = 2) +
    coord_equal(),
  ggraph(lo) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size),
                    data = function(df) {
                      filter(df, vertex %in% david.v2p[.("A", "US (post)")]$vertex)
                    },
                    fill = "purple", shape = 21) +
    scale_size_area(max_size = 2) +
    labs(title = "A, post-travel") +
    coord_equal(),
  ggraph(lo) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size),
                    data = function(df) {
                      filter(df, vertex %in% david.v2p[.("B", "pre-Salmonella")]$vertex)
                    },
                    fill = "red", shape = 21) +
    scale_size_area(max_size = 2) +
    labs(title = "B, pre-Salmonella") +
    coord_equal(),
  ggraph(lo) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size),
                    data = function(df) {
                      filter(df, vertex %in% david.v2p[.("B", "post-Salmonella")]$vertex)
                    },
                    fill = "orange", shape = 21) +
    scale_size_area(max_size = 2) +
    labs(title = "B, post-Salmonella") +
    coord_equal()
)
events <- last_plot()


# basins ------------------------------------------------------------------


ggraph(lo) +
  geom_edge_link0() +
  geom_node_point(aes(fill = as.factor(basin), size = size),
                  data = function(df) filter(df, !in.singleton), shape = 21) +
  labs(fill = "basin") +
  guides(size = FALSE) +
  theme_graph(base_family = "Helvetica") +
  coord_equal()
basins <- last_plot()


# basin series -----------------------------------------------


david.v2p[, basin := as.factor(basin)]
david.sample.basins <- david.v2p[, .N, by = .(subject, day, point, event, basin)]
theme_set(theme_cowplot())
david.sample.basins %>%
  ggplot(aes(x = day, y = basin)) +
  geom_tile(data = function(dt) dt[is.na(basin)], fill = "grey50") +
  geom_tile(aes(fill = event), data = function(dt) dt[!is.na(basin)]) +
  scale_fill_brewer(palette = "Set1") +  
  scale_y_discrete(drop = FALSE) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(subject ~ .)
series <- last_plot()

# distribs ----------------------------------------------------------------


david.sample.basins %>%
  group_by(subject, event) %>%
  mutate(frac = N / sum(N)) %>%
  group_by(subject, event, basin) %>%
  summarize(frac = sum(frac)) %>%
  ggplot(aes(x = event)) +
  geom_col(aes(y = frac, group = basin, fill = basin),
           position = "stack") +
  coord_flip() +
  labs(y = "fraction samples") +
  facet_grid(subject ~ ., scales = "free")
distribs <- last_plot()

# correlation function ----------------------------------------------------


david.persistence <- david.v2p %>%
  split(by = c("subject", "event")) %>%
  lapply(function(df) basin.persistence(df$basin, df$day, scale = TRUE)) %>%
  rbindlist(idcol = "subject.event") %>%
  .[, c("subject", "event") := tstrsplit(subject.event, "\\.")] %>%
  .[, subject.event := NULL] %>%
  .[, basin := as.factor(as.numeric(basin))] %>%
  # .[, N := .N, by = .(basin, delta.t, subject, event)]
  .[, N := uniqueN(delta.t), by = .(basin, subject, event)]
setkey(david.persistence, subject)
pl <- lapply(split(david.persistence[!is.na(basin)], by = "subject"),
          function(dt) {
            dt %>%
              ggplot(aes(x = delta.t, y = f, color = basin)) +
              geom_point(data = function(df) filter(df, N <= 10),
                          alpha = 0.5) +
              geom_smooth(aes(color = basin, fill = basin),
                          data = function(df) filter(df, N > 10)) +
              scale_color_hue(drop = FALSE) +
              scale_fill_hue(drop = FALSE) +
              coord_cartesian(ylim = c(0, 1)) +
              scale_y_continuous(breaks = c(0, 0.5, 1)) +
              labs(x = "interval (days)", y = "correlation") +
              facet_wrap(~ event, ncol = 1)
            })
theme_set(theme_cowplot() + theme(legend.position = "none",
                                  axis.text = element_text(size = title.size),
                                  axis.title = element_text(size = title.size),
                                  strip.text = element_text(size = title.size)))
plot_grid(pl[[1]], pl[[2]], nrow = 1)
correlation <- last_plot()


# compiled ----------------------------------------------------------------

plot_grid(fsubject, events, ncol = 2, labels = "AUTO")
save_plot("../../figures/tda/paper/fig3.pdf", last_plot(), ncol = 2,
          base_width = 4)
plot_grid(basins + guides(fill = guide_legend(nrow = 10)),
          plot_grid(distribs + theme(legend.position = "none",
                                     axis.text = element_text(size = title.size),
                                     axis.title = element_text(size = title.size)), 
                    correlation, 
                    ncol = 2, labels = c("B", "C")),
          series + theme_cowplot(font_size = title.size) + 
            guides(fill = guide_legend(direction = "horizontal")) +
            theme(legend.position = "bottom", axis.text.y = element_blank(),
                  axis.ticks = element_blank()),
          labels = c("A", "", "D"),
          ncol = 1,
          rel_heights = c(2, 2, 1)
          )
save_plot("../../figures/tda/paper/fig4.pdf", last_plot(), nrow = 3, base_width = 8)
