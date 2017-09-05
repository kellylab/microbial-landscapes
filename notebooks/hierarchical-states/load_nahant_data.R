library(tidyverse)
library(data.table)

source("~/github_bitbucket/wkc.r-utils/df_2_mat.R")


data.dir <- "../../data/"
read_nahant <- function(fn, data.dir) {
  # TODO complete this
  dat <- fread(paste0(data.dir, fn), header = TRUE)
  dat <- melt(dat, id.vars = c("OTU", "ConsensusLineage"),
                     variable.name = "day")
  dat[, day := as.numeric(as.character(day))]
  dat[, value := as.numeric(value)]
  # merge OTU counts by consensus lineage
  dat[, .(freq = sum(value)), by = .(ConsensusLineage, day)]
  # parse lineage
  # nothing is mapped to species resolution
  dat[, ConsensusLineage := gsub(";", "", ConsensusLineage)]
  dat[, c("foo", "kingdom", "phylum", "class", "order", "family",
                 "genus") := tstrsplit(ConsensusLineage,
                                       "[[:alpha:]]__")]
  dat[, foo := NULL]
}
nahant <- lapply(list(pro = "BactNorm.txt", euk = "EukNorm.txt"), read_nahant,
                 data.dir = data.dir) %>%
  rbindlist
