library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(cowplot)

scripts.dir <- "../r/"
source(paste0(scripts.dir, "load-prochlorococcus-data.R"))
source("k.first.R")
source("vertex.2.points.R")
source("mapper.2.igraph.R")
source("plot.mapper.R")

prochlorococcus[, sample := paste(site, cruiseid, depth, sep = "-")]
prochlorococcus[, total.abundance := sum(abundance), by = sample]
prochlorococcus <- prochlorococcus[total.abundance > 0]
prochlorococcus[, freq := abundance / total.abundance, by = sample]
distribs <- dcast(prochlorococcus, sample ~ ecotype, value.var = "freq")
samples <- distribs[, 1]
jsds <- JSD(as.matrix(distribs[, -1]))
