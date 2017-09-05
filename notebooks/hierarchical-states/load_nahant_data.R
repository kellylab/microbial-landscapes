library(tidyverse)
library(data.table)

scripts.dir <- "../../scripts/r/"
source(paste0(scripts.dir, "lineage_2_taxonomy.R"))
# make this a real package someday
source("~/github_bitbucket/wkc.r-utils/df_2_mat.R")

data.dir <- "../../data/"
read_nahant <- function(fn, data.dir) {
  dat <- fread(paste0(data.dir, fn), header = TRUE)
  # lineage columns are named differently
  if (fn == "BactNorm.txt") {
    setnames(dat, "ConsensusLineage", "Lineage")
  }
  dat <- melt(dat, id.vars = c("OTU", "Lineage"),
                     variable.name = "day")
  dat[, day := as.numeric(as.character(day))]
  dat[, value := as.numeric(value)]
  # merge OTU counts by consensus lineage
  dat[, .(freq = sum(value)), by = .(Lineage, day)]
  # parse lineage
  dat[, c("kingdom", 
          "phylum", 
          "class", 
          "order", 
          "family", 
          "genus", 
          "species") := lineage_2_taxonomy(Lineage), by = Lineage]
}
nahant <- lapply(list(pro = "BactNorm.txt", euk = "EukNorm.txt"), read_nahant,
                 data.dir = data.dir) %>%
  rbindlist
