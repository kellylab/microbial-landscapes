library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(cowplot)
library(combinat)

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
samples <- distribs[, sample]
# jsds <- JSD(as.matrix(distribs[, -1])) # currently not work
distribs <- as.matrix(distribs[, -1])
rownames(distribs) <- samples
jsds <- data.table(t(combn(samples, 2)))
names(jsds) <- c("sample.i", "sample.j")
jsds[, jsd := JSD(distribs[c(sample.i, sample.j),], unit = "log2"),
     by = 1:nrow(jsds)]
