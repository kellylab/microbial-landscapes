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
  geom_jitter(size = 0.1) +
  geom_smooth(aes(fill = basin)) +
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
ggplot(david.persistence, aes(x = delta.t, y = f, color = basin)) +
  geom_smooth(aes(fill = basin)) +
  facet_wrap(~ subject, nrow = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_hue(drop = FALSE) +
  scale_fill_hue(drop = FALSE) +
  labs(x = "interval (days)", y = "correlation")
