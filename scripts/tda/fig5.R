#' # Basin persistence time
# setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(cowplot)

scripts.dir <- "../r/"
mapper.dir <- "../../output/mapper/"
figs.dir <- "../../figures/tda/paper/"
source("basin.persistence.R")

read.mapper.graph <- function(directory) {
  vertices <- fread(paste0(directory, "/vertices.txt"))
  edges <- fread(paste0(directory, "/edges.txt"))
  graf <- tbl_graph(vertices, edges, FALSE)
  v2p <- fread(paste0(directory, "/vertices-to-points.txt"))
  list(graph = graf, mapping = v2p)
}
merge.mpr.samples <- function(mpr, dt, by.y = "sample") {
  vertices <- mpr$graph %>%
    activate(nodes) %>%
    as.data.table %>%
    merge(mpr$mapping, by = "vertex") %>%
    merge(dt, by.x = "point", by.y = by.y)
  vertices
}


# fig 1 -------------------------------------------------------------------
source(paste0(scripts.dir, "load_cholera_data.R"))
cholera.mapper <- read.mapper.graph(paste0(output.dir, "cholera"))
cholera.v2p <- merge.mpr.samples(cholera.mapper, gordon.samples)

subplots <- list() # list of subplots

set.seed(1)
lo <- create_layout(cholera.mapper$graph, "fr", niter = 1000)
xy <- lo[, c("x", "y")]
# color by fraction samples marked diarrhea
theme_set(theme_graph(base_family = "Helvetica"))
ggraph(cholera.mapper$graph, "manual", 
                          node.positions = xy) +
  geom_edge_link0() +
  geom_node_point(aes(color = f.state, size = size),
                  data = function(df) filter(df, in.singleton),
                  fill = "white", shape = 21) +
  geom_node_point(aes(fill = f.state, size = size),
                  data = function(df) filter(df, !in.singleton), shape = 21) +
  scale_color_distiller(palette = "Spectral") +
  scale_fill_distiller(palette = "Spectral") +
  labs(color = "fraction\ndiarrhea") +
  guides(fill = FALSE) +
  theme_graph(base_family = "Helvetica") +
  coord_equal()
subplots$fstate <- last_plot()
ggraph(cholera.mapper$graph, "manual", node.positions = xy) +
  geom_edge_link0() +
  geom_node_point(aes(size = size, color = scaled.knn)) +
  scale_color_distiller(palette = "Greys") +
  labs(color = "inverse\ndensity") +
  theme_graph(base_family = "Helvetica") +
  coord_equal() +
  guides(size = FALSE) +
  scale_edge_color_distiller(palette = "Greys")
subplots$density <- last_plot()
ggraph(cholera.mapper$graph, "manual", node.positions = xy) +
  geom_edge_link0() +
  geom_node_point(aes(fill = as.factor(basin), size = size),
                  # data = function(df) filter(df, !in.singleton), 
                  shape = 21) +
  labs(fill = "basin") +
  scale_fill_brewer(palette = "Accent") +
  theme_graph(base_family = "Helvetica") +
  coord_equal() +
  guides(size = FALSE) 
subplots$basin <- last_plot()

#' # Subject dynamics

#' ## Subject dynamics in terms of basins
theme_set(theme_cowplot())
basins <- as.data.frame(activate(cholera.mapper$graph, nodes))$basin %>%
  unique %>%
  as.numeric %>%
  sort(na.last = TRUE) %>%
  as.character
p2basin <- cholera.mapper$graph %>%
  activate(nodes) %>%
  as.data.frame %>%
  select(vertex, basin) %>%
  merge(v2p, by = "vertex") %>%
  as.data.table %>%
  .[, .N, by = .(subject, diagnosis, id, hour, basin)]
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
  theme(axis.ticks.y = element_blank()) +
  labs(x = "hour", fill = "basin") +
  scale_fill_brewer(palette = "Accent", drop = FALSE) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), strip.placement = "outside")
d2basin <- p2basin[, .(N = sum(N)), by = .(subject, diagnosis, basin)]
d2basin[, frac := N / sum(N), by = .(subject, diagnosis)]
setkey(d2basin, diagnosis)
d2basin["diarrhea"] %>%
  ggplot(aes(x = subject, y = frac)) +
  geom_col(aes(fill = basin)) +
  coord_flip() +
  scale_fill_brewer(palette = "Accent", drop = FALSE) +
  facet_grid(subject ~ ., scales = "free_y") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        strip.background = element_blank(), strip.text = element_blank()) +
  labs(fill = "basin", y = "fraction samples") #+
subplots$basin.distribs <- last_plot()
# plot_grid(subplots$fstate, subplots$knn, subplots$basins,
#           plot_grid(subplots$basin.series, subplots$basin.distribs, ncol = 2,
#                     rel_widths = c(1, 1), align = "hv", axis = "lt"),
#           labels = "AUTO",
#           hjust = c(-0.5, -0.5, -0.5, 0.5),
#           vjust = c(1.5, 1.5, 1.5, 0))
# save_plot(paste0(figs.dir, "paper/fig2.pdf"),
#           last_plot(),
#           base_aspect_ratio = 1.3, base_width = 8,
#           ncol = 1, nrow = 2)


# fig 2 -------------------------------------------------------------------

source(paste0(scripts.dir, "load_david_data.R"))
david.mapper <- read.mapper.graph(paste0(mapper.dir, "david/"))
david.v2p <- merge.mpr.samples(david.mapper, david.samples)
# the line below is a kludge
setnames(david.v2p, c("subject.x", "subject.y"), c("f.subject", "subject"))
set.seed(0)
david.xy <- create_layout(david.mapper$graph, "fr", niter = 1000)[, c("x", "y")]
david.subplots <- list()
ggraph(david.mapper$graph, "manual", node.positions = david.xy) +
  geom_edge_link0() +
  geom_node_point(aes(color = f.subject, size = size),
                  data = function(df) filter(df, in.singleton),
                  fill = "white", shape = 21) +
  geom_node_point(aes(fill = subject, size = size),
                  data = function(df) filter(df, !in.singleton), shape = 21) +
  labs(color = "fraction A") +
  guides(fill = FALSE) +
  scale_color_distiller(palette = "Spectral") +
  scale_fill_distiller(palette = "Spectral") +
  theme_graph(base_family = "Helvetica") +
  coord_equal()
david.subplots$fsubject <- last_plot()

# events
setkey(david.v2p, subject, event)
plot_grid(
  ggraph(david.mapper$graph, "manual", node.positions = david.xy) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size), 
                    data = function(df) {
                      filter(df, vertex %in% david.v2p[.("A", "US (pre)")]$vertex)
                    },
                    fill = "blue", shape = 21) +
    labs(title = "A, pre-travel") +
    theme_graph(base_family = "Helvetica") +
    coord_equal(),
  ggraph(david.mapper$graph, "manual", node.positions = david.xy) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size), 
                    data = function(df) {
                      filter(df, vertex %in% david.v2p[.("A", "US (post)")]$vertex)
                    },
                    fill = "purple", shape = 21) +
    labs(title = "A, post-travel") +
    theme_graph(base_family = "Helvetica") +
    coord_equal(),
  ggraph(david.mapper$graph, "manual", node.positions = david.xy) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size), 
                    data = function(df) {
                      filter(df, vertex %in% david.v2p[.("B", "pre-Salmonella")]$vertex)
                    },
                    fill = "red", shape = 21) +
    labs(title = "B, pre-Salmonella") +
    theme_graph(base_family = "Helvetica") +
    coord_equal(),
  ggraph(david.mapper$graph, "manual", node.positions = david.xy) +
    geom_edge_link0() +
    geom_node_point(aes(size = size), color = "grey") +
    geom_node_point(aes(size = size), 
                    data = function(df) {
                      filter(df, vertex %in% david.v2p[.("B", "post-Salmonella")]$vertex)
                    },
                    fill = "orange", shape = 21) +
    labs(title = "B, post-Salmonella") +
    theme_graph(base_family = "Helvetica") +
    coord_equal()
)


ggraph(david.mapper$graph, "manual", node.positions = david.xy) +
  geom_edge_link0() +
  geom_node_point(aes(size = size),
                  data = function(df) filter(df, in.singleton),
                  fill = "white", shape = 21) +
  geom_node_point(aes(fill = as.factor(basin), size = size),
                  data = function(df) filter(df, !in.singleton), shape = 21) +
  labs(fill = "basin") +
  guides(size = FALSE) +
  # geom_node_point(color = "black",
  #                 data = function(df) filter(df, is.extremum)) +
  theme_graph(base_family = "Helvetica") +
  coord_equal()
david.subplots$basin <- last_plot()

david.v2p[, basin := as.factor(basin)]
david.sample.basins <- david.v2p[, .N, by = .(subject, day, point, event, basin)] 
david.sample.basins %>%
  ggplot(aes(x = day, y = basin)) +
  geom_tile(data = function(dt) filter(dt, is.na(basin)), fill = "grey50") +
  geom_tile(aes(fill = event), data = function(dt) filter(dt, !is.na(basin))) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(subject ~ .)
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
# fig 3 --------------------------------------------------------------------

source(paste0(scripts.dir, "load-prochlorococcus-data.R"))
hotbats.mapper <- read.mapper.graph(paste0(mapper.dir, "prochlorococcus"))
hotbats.v2p <- merge.mpr.samples(hotbats.mapper, prochlorococcus.samples)
set.seed(1)
hotbats.xy <- create_layout(hotbats.mapper$graph, "fr", niter = 1000)
ggraph(hotbats.mapper$graph, "manual", node.positions = hotbats.xy) +
  geom_edge_link0() +
  geom_node_point(aes(size = size, fill = mean.depth), shape = 21) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  guides(fill = guide_colorbar(reverse = TRUE)) +
  theme_graph(base_family = "Helvetica") +
  labs(fill = "depth (m)") +
  coord_equal()
ggraph(hotbats.mapper$graph, "manual", node.positions = hotbats.xy) +
  geom_edge_link0() +
  geom_node_point(aes(size = size, fill = mean.depth), shape = 21) +
  scale_fill_distiller(palette = "Spectral") +
  theme_graph(base_family = "Helvetica") +
  labs(fill = "temp (C)") +
  coord_equal()
hotbats.mapper$graph <- hotbats.mapper$graph %>% 
  activate(nodes) %>% 
  mutate(basin = as.factor(basin)) %>% 
  mutate(basin = mapply(function(b, x) if (x) NA else b), b = basin, x = in.singleton)
ggraph(hotbats.mapper$graph, "manual", node.positions = hotbats.xy) +
  geom_edge_link0() +
  geom_node_point(aes(size = size, fill = basin), shape = 21) +
  theme_graph(base_family = "Helvetica") +
  coord_equal()




# cholera -----------------------------------------------------------------

# source(paste0(scripts.dir, "load_cholera_data.R"))
# cholera.mpr <- read.mapper.graph(paste0(mapper.dir, "cholera/"))
# cholera.v2p <- merge.mpr.samples(cholera.mpr, gordon.samples)
cholera.persistence <- split(cholera.v2p, by = c("subject", "diagnosis")) %>%
  lapply(function(df) basin.persistence(df$basin, df$hour, scale = TRUE)) %>%
  rbindlist(idcol = "group") %>%
  .[, c("subject", "diagnosis") := tstrsplit(group, "\\.")] %>%
  .[, basin := as.factor(as.numeric(basin))]
setkey(cholera.persistence, diagnosis)
ggplot(cholera.persistence["diarrhea"],
                   aes(x = delta.t, y = f, color = basin)) +
  geom_smooth(aes(fill = basin)) +
  stat_summary(geom = "point", fun.y = mean, size = 0.1) +
  scale_color_brewer(palette = "Accent", drop = FALSE) +
  scale_fill_brewer(palette = "Accent", drop = FALSE) +
  facet_wrap(~ subject) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "interval (hours)", y = "correlation")
pcholera <- last_plot()


# david -------------------------------------------------------------------


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
              labs(x = "interval (days)", y = "correlation") +
              facet_wrap(~ event, ncol = 1) 
            })
plot_grid(pl[[1]] + theme(legend.position = "none"),
          pl[[2]] + theme(legend.position = "none"),
          get_legend(pl[[1]]),
          nrow = 1, rel_widths = c(4, 4, 1), labels = c("A", "B", "")
          )
pdavid <- ggplot(david.persistence[!is.na(basin)],
                 aes(x = delta.t, y = f, color = basin)) +
  geom_smooth(aes(fill = basin)) +
  stat_summary(geom = "point", fun.y = mean, size = 0.1) +
  facet_wrap(~ subject, ncol = 2) +
  # coord_cartesian(ylim = c(0, 1)) +
  scale_color_hue(drop = FALSE) +
  scale_fill_hue(drop = FALSE) +
  labs(x = "interval (days)", y = "correlation")


# prochlorococcus --------------------------------------------------------

hotbats.persistence <- hotbats.v2p %>%
  split(by = c("site", "depth")) %>%
  lapply(function(df) basin.persistence(df$basin, df$month, scale = TRUE)) %>%
  rbindlist(idcol = "site.depth") %>%
  .[, c("site", "depth") := tstrsplit(site.depth, "\\.")] %>%
  .[, site.depth := NULL] %>%
  .[, basin := as.factor(as.numeric(basin))] %>%
  .[, depth := as.factor(as.numeric(depth))]
pl <- hotbats.persistence[!is.na(basin)] %>%
  split(by = "site") %>%
  lapply(function(df) {
    ggplot(df, aes(x = delta.t, y = f, color = basin)) +
      stat_smooth(aes(fill = basin), data = function(dt) {
        # filter out depth-basin combos that have < 10 unique delta.t values
        # these mess up the smoothing calculation for the entire facet!
        dt %>%
          group_by(depth, basin) %>%
          mutate(nint = uniqueN(delta.t)) %>%
          ungroup %>%
          filter(nint >= 10)
      }) +
      stat_summary(geom = "point", fun.y = mean, size = 0.1) +
      facet_wrap(~ depth) +
      scale_color_hue(drop = FALSE) +
      scale_fill_hue(drop = FALSE) +
      labs(x = "interval (months)", y = "correlation")
  })
plot_grid(plot_grid(pl[[1]] + theme(legend.position = "none"),
                    pl[[2]] + theme(legend.position = "none"),
                    nrow = 2, labels = toupper(names(pl))),
  get_legend(pl[[1]]),
  rel_widths = c(8, 1)
)
pchloro <- last_plot()

# put it all together -----------------------------------------------------

plot_grid(plot_grid(pcholera, pdavid, nrow = 1, labels = "AUTO"),
          pchloro, nrow = 2, labels = c("", "C"))
save_plot(paste0(figs.dir, "fig5.pdf"), last_plot(), nrow = 2, ncol = 1,
          base_width = 8)
