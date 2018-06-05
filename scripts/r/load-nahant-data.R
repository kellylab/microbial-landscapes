library(tidyverse)
library(data.table)

data.dir <- "../../data/"

lineage.2.taxonomy <- function(lineage) {
  llist <- data.table::tstrsplit(lineage, ";")
  lapply(llist, function(v) {
    sapply(v, function(w) {
      out <- sub("[[:alpha:]]__", "", w)
      if (out == "") {
        out <- NA_character_
      }
      out
    })
  })
}
read.nahant <- function(fn, data.dir) {
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
  #dat[, .(freq = sum(value)), by = .(Lineage, day)]
  # parse lineage
  dat[, c("kingdom",
          "phylum",
          "class",
          "order",
          "family",
          "genus",
          "species") := lineage.2.taxonomy(Lineage), by = Lineage]
}
nahant <- lapply(list(pro = "BactNorm.txt", euk = "EukNorm.txt"), read.nahant,
                 data.dir = data.dir) %>%
  rbindlist
