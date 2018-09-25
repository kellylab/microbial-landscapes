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

# cholera -----------------------------------------------------------------

source(paste0(scripts.dir, "load_cholera_data.R"))
cholera.mpr <- read.mapper.graph(paste0(mapper.dir, "cholera/"))
cholera.v2p <- merge.mpr.samples(cholera.mpr, gordon.samples)
cholera.persistence <- split(cholera.v2p, by = c("subject", "diagnosis")) %>%
  lapply(function(df) basin.persistence(df$basin, df$hour, scale = TRUE)) %>%
  rbindlist(idcol = "group") %>%
  .[, c("subject", "diagnosis") := tstrsplit(group, "\\.")] %>%
  .[, basin := as.factor(as.numeric(basin))]
setkey(cholera.persistence, diagnosis)
pcholera <- ggplot(cholera.persistence["diarrhea"],
                   aes(x = delta.t, y = f, color = basin)) +
  geom_smooth(aes(fill = basin)) +
  stat_summary(geom = "point", fun.y = mean, size = 0.1) +
  scale_color_hue(drop = FALSE) +
  scale_fill_hue(drop = FALSE) +
  facet_wrap(~ subject) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "interval (hours)", y = "correlation")


# david -------------------------------------------------------------------


source(paste0(scripts.dir, "load_david_data.R"))
david.mpr <- read.mapper.graph(paste0(mapper.dir, "david/"))
david.v2p <- merge.mpr.samples(david.mpr, david.samples)
setnames(david.v2p, c("subject.x", "subject.y"), c("frac.a", "subject")) #KLUDGE
david.persistence <- david.v2p %>%
  split(by = "subject") %>%
  lapply(function(df) basin.persistence(df$basin, df$day, scale = TRUE)) %>%
  rbindlist(idcol = "subject") %>%
  .[, basin := as.factor(as.numeric(basin))]
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

source(paste0(scripts.dir, "load-prochlorococcus-data.R"))
hotbats.mpr <- read.mapper.graph(paste0(mapper.dir, "prochlorococcus/"))
hotbats.v2p <- merge.mpr.samples(hotbats.mpr, prochlorococcus.samples)
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
pchloro <- plot_grid(plot_grid(pl[[1]] + theme(legend.position = "none"),
                    pl[[2]] + theme(legend.position = "none"),
                    nrow = 2, labels = toupper(names(pl))),
  get_legend(pl[[1]]),
  rel_widths = c(8, 1)
)

# put it all together -----------------------------------------------------

plot_grid(plot_grid(pcholera, pdavid, nrow = 1, labels = "AUTO"),
          pchloro, nrow = 2, labels = c("", "C"))
save_plot(paste0(figs.dir, "fig5.pdf"), last_plot(), nrow = 2, ncol = 1,
          base_width = 8)
