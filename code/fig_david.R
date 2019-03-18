# setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(cowplot)

mapper.dir <- "mapper-output/david/"
utils.dir <- "utils/"
figs.dir <- "../figures/"
for (script in list.files(utils.dir, full.names = TRUE)) source(script)
source("load_david_data.R")
david.mapper <- read.mapper.graph(mapper.dir)
david.v2p <- merge.mpr.samples(david.mapper, david.samples)
plotter <- function(v) {
  plot.mapper.graph(david.mapper$graph,
                    node = geom_node_point(aes_(size = ~size, fill = v),
                                           shape = 21),
                    seed = 13,
                    exclude.singletons = TRUE) +
    guides(size = FALSE) +
    scale_size_area(max_size = 4)
}

# frac subject ------------------------------------------------------------

plotter(~f.subject) +
  scale_fill_distiller(palette = "Spectral", breaks = c(0, 0.5, 1),
                       labels = c("all B", "", "all A")) +
  labs(fill = "subject")
fsubject <- last_plot()


# events ------------------------------------------------------------------


set.seed(13)
lo <- david.mapper$graph %>% # mapper graph with no outlier
  activate(nodes) %>%
  filter(!in.singleton) %>%
  create_layout("fr", niter = 1000)
title.size <- 10
base.size <- 8
plt.event <- function(subj, ev, color, label) {
  setkey(david.v2p, subject, event)
  vs <- david.v2p[.(subj, ev)]$vertex
  ggraph(lo) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size),
                    data = function(df) filter(df, vertex %in% vs),
                    fill = color, shape = 21) +
    labs(title = label) +
    theme_graph(base_family = "Helvetica", title_size = title.size,
                      base_size = base.size) +
    theme(legend.position = "none", plot.margin = margin(1, 1, 1, 1)) +
    scale_size_area(max_size = 2) +
    coord_equal()
}
plot_grid(plotlist = mapply(plt.event,
                            subj = c("A", "A", "B", "B"),
                            ev = c("US (pre)", "US (post)",
                                      "pre-Salmonella", "post-Salmonella"
                                      ),
                            color = c("blue", "purple", "red", "orange"),
                            label = c("A, pre-travel", "A, post-travel",
                                      "B, pre-Salmonella", "B, post-Salmonella"
                                      ),
                            SIMPLIFY = FALSE)
          )
events <- last_plot()


# basins ------------------------------------------------------------------

plotter(~as.factor(basin)) + labs(fill = "basin")
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
  split(by = c("subject", "event", "healthy")) %>%
  lapply(function(df) basin.persistence(df$basin, df$day, scale = TRUE)) %>%
  rbindlist(idcol = "id") %>%
  .[, basin := as.factor(as.numeric(basin))] %>%
  .[, N := uniqueN(delta.t), by = .(id, basin)] %>%
  .[, c("subject", "event", "healthy") := tstrsplit(id, "\\.",
                                                    type.convert = TRUE)] %>%
  .[, id := NULL]
theme_set(theme_cowplot() + theme(legend.position = "none",
                                  axis.text = element_text(size = title.size),
                                  axis.title = element_text(size = title.size),
                                  strip.text = element_text(size = title.size)))
graphic.size <- 0.5
ggplot(david.persistence[!is.na(basin) & healthy],
       aes(x = delta.t, y = f, color = basin)) +
  geom_point(data = function(df) filter(df, N <= 10),
             alpha = 0.5, size = graphic.size) +
  geom_smooth(aes(color = basin, fill = basin),
              data = function(df) filter(df, N > 10),
              size = graphic.size) +
  scale_color_hue(drop = FALSE) +
  scale_fill_hue(drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "interval (days)", y = "correlation") +
  facet_grid(paste(subject, event, sep = ", ") ~ .) +
  theme(strip.text.y = element_text(angle = 0))
correlation <- last_plot()


# compiled ----------------------------------------------------------------

plot_grid(fsubject, events,
          basins, plot_grid(distribs, correlation, nrow = 2),
          ncol = 2, rel_heights = c(2, 3), labels = "AUTO")
save_plot(paste0(figs.dir, "/fig3.pdf"), last_plot(), ncol = 2,
          base_width = 4, base_height = 8)
