
# setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(cowplot)

utils.dir <- "utils/"
scripts.dir <- "../r/"
output.dir <- "../../output/mapper/"
for (script in list.files(utils.dir, full.names = TRUE)) source(script)
source(paste0(scripts.dir, "load_cholera_data.R"))
cholera.mapper <- read.mapper.graph(paste0(output.dir, "cholera"))
cholera.v2p <- merge.mpr.samples(cholera.mapper, gordon.samples)
# lay out Mapper graph without singletons
set.seed(1)
lo <- cholera.mapper$graph %>%
  activate(nodes) %>%
  filter(!in.singleton) %>%
  create_layout("fr", niter = 1000)
theme_set(theme_graph(base_family = "Helvetica"))
# color settings for basins of attraction
basin.palette <- "Set3"
na.color <- "grey50"
subplots <- list() # list of subplots

# fraction diarrhea -------------------------------------------------------


ggraph(lo) + # color by fraction samples marked diarrhea
  geom_edge_link0() +
  geom_node_point(aes(fill = f.state, size = size),
                  data = function(df) filter(df, !in.singleton), shape = 21) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill = "fraction\ndiarrhea") +
  guides(size = FALSE) +
  theme_graph(base_family = "Helvetica") +
  coord_equal()
subplots$fstate <- last_plot()

# density -----------------------------------------------------------------


ggraph(lo) +
  geom_edge_link0() +
  geom_node_point(aes(size = size, color = scaled.knn)) +
  scale_color_distiller(palette = "Greys") +
  labs(color = "Q") +
  theme_graph(base_family = "Helvetica") +
  coord_equal() +
  guides(size = FALSE) +
  scale_edge_color_distiller(palette = "Greys")
subplots$density <- last_plot()


# basins ------------------------------------------------------------------


ggraph(lo) +
  geom_edge_link0() +
  geom_node_point(aes(fill = as.factor(basin), size = size),
                  # data = function(df) filter(df, !in.singleton),
                  shape = 21) +
  labs(fill = "basin") +
  scale_fill_brewer(palette = basin.palette, na.value = na.color) +
  theme_graph(base_family = "Helvetica") +
  coord_equal() +
  guides(size = FALSE)
subplots$basin <- last_plot()


# basin time series -------------------------------------------------------


theme_set(theme_cowplot())
basins <- as.data.frame(activate(cholera.mapper$graph, nodes))$basin %>%
  unique %>%
  as.numeric %>%
  sort(na.last = TRUE) %>%
  as.character
p2basin <- cholera.v2p[, .N, by = .(subject, diagnosis, id, hour, basin)]
p2basin[, time := mapply(function(diagnosis, id, hour) {
  if (diagnosis == "diarrhea") {
    hour
  } else {
    as.numeric(strsplit(id, "d")[[1]][[2]])
  }
}, diagnosis = diagnosis, id = id, hour = hour)]
p2basin[, time.unit := sapply(diagnosis, function(x) {
  if (x == "diarrhea") "hour" else "day"
})]
p2basin[, i := rank(basin), by = .(subject, id)]
p2basin[, basin := as.factor(basin)]
theme_set(theme_cowplot(font_size = 10))
setkey(p2basin, diagnosis)
p2basin["diarrhea"] %>%
  ggplot(aes(x = time, y = basin)) +
  geom_tile(aes(fill = basin)) +
  facet_grid(subject ~ ., scales = "free_y", switch = "y") +
  theme(axis.ticks.y = element_blank(), panel.spacing = unit(1, "points")) +
  labs(x = "hour", fill = "basin") +
  scale_fill_brewer(palette = basin.palette, na.value = na.color,
                    drop = FALSE) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), strip.placement = "outside")
subplots$basin.series <- last_plot()

# basin distribution ------------------------------------------------------


d2basin <- p2basin[, .(N = sum(N)), by = .(subject, diagnosis, basin)]
d2basin[, frac := N / sum(N), by = .(subject, diagnosis)]
setkey(d2basin, diagnosis)
d2basin["diarrhea"] %>%
  ggplot(aes(x = subject, y = frac)) +
  geom_col(aes(fill = basin)) +
  coord_flip() +
  scale_fill_brewer(palette = basin.palette, na.value = na.color, drop = FALSE) +
  facet_grid(subject ~ ., scales = "free_y") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), panel.spacing = unit(1, "points")) +
  labs(fill = "basin", y = "fraction samples") #+
subplots$basin.distribs <- last_plot()
subplots$basin.combined <- plot_grid(subplots$basin.series +
                                       theme(legend.position = "none"),
                                     subplots$basin.distribs +
                                       theme(strip.background = element_blank(),
                                             strip.text = element_blank(),
                                             legend.position = "none"),
                                     align = "hv", axis = "t"
                                     )

# basin correlation function ----------------------------------------------


cholera.persistence <- split(cholera.v2p, by = c("subject", "diagnosis")) %>%
  lapply(function(df) basin.persistence(df$basin, df$hour, scale = TRUE)) %>%
  rbindlist(idcol = "group") %>%
  .[, c("subject", "diagnosis") := tstrsplit(group, "\\.")] %>%
  .[, basin := as.factor(as.numeric(basin))]
setkey(cholera.persistence, diagnosis)
cholera.persistence["diarrhea"] %>%
  .[!is.na(basin)] %>%
  ggplot(aes(x = delta.t, y = f, color = basin)) +
  geom_smooth(aes(fill = basin)) +
  stat_summary(geom = "point", fun.y = mean, size = 0.1) +
  scale_color_brewer(palette = basin.palette, na.value = na.color, drop = FALSE) +
  scale_fill_brewer(palette = basin.palette, na.value = na.color,  drop = FALSE) +
  facet_wrap(~ subject, nrow = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "interval (hours)", y = "correlation")
subplots$correlation <- last_plot()


# combined figure ---------------------------------------------------------


plot_grid(plot_grid(subplots$fstate, subplots$basin,
                    labels = "AUTO", align = "h", nrow = 1),
          subplots$basin.combined,
          subplots$correlation + theme(legend.position = "none"),
          rel_heights = c(2, 1, 1),
          labels = c("", "C", "D"),
          label_y = c(1, 1.2, 1),
          ncol = 1)
save_plot(last_plot(), filename = "../../figures/tda/paper/fig2.pdf",
          ncol = 1, base_width = 8, base_height = 8)
