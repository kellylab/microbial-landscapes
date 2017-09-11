library(tidyverse)
library(data.table)
library(philentropy)
#library(infograf)

source("../../../../wkc.r-utils/df_2_mat.R")
source("../../../../wkc.r-utils/mat_2_df.R")

# general function to get JSDs from pairwise data
jsds_from_d <- function(d, sample.cn, var.cn, val.cn, suf = c(".i", ".j")) {
  form <- paste(sample.cn, "~", var.cn)
  d <- dcast(d, form, value.var = val.cn, fill = 0)
  rn <- d[[sample.cn]]
  d <- select(d, -matches(sample.cn))
  m <- as.matrix(d)
  rownames(m) <- rn
  jsds <- JSD(m)
  rownames(jsds) <- rn
  colnames(jsds) <- rn
  si <- paste0(sample.cn, suf[1])
  sj <- paste0(sample.cn, suf[2])
  jsds <- mat_2_df(jsds, row.names = si)
  jsds <- as.data.table(jsds)
  melt(jsds,
       id.vars       = si,
       variable.name = sj,
       value.name = "jsd",
       variable.factor = FALSE)
}

# Prochlorococcus ---------------------------------------------------------

# load data
source("load_prochlorococcus_data.R")
# for multiple samples in a month take mean abundance
prochlorococcus <- prochlorococcus[, .(abundance = mean(abundance)),
                                   by = .(site,
                                          year,
                                          month,
                                          cal.month,
                                          depth,
                                          ecotype)]
# relative abundance
prochlorococcus[, rel.abund := abundance / sum(abundance),
                by = .(site, month, depth)]
# get rid of rows where relative abundance is NaN; this is when nothing was
# sampled
prochlorococcus <- prochlorococcus[!is.nan(rel.abund)]
prochlorococcus[, sample := paste(year, cal.month, depth, sep = "_")]
# pairwise JSDs
proch.jsds <- prochlorococcus[, jsds_from_d(.SD,
                                            "sample",
                                            "ecotype",
                                            "rel.abund"),
                              by = site]
proch.jsds[, c("year.i", "month.i", "depth.i") := tstrsplit(sample.i, "_")]
proch.jsds[, c("year.j", "month.j", "depth.j") := tstrsplit(sample.j, "_")]
proch.jsds[, c("sample.i", "sample.j") := NULL]
fwrite(proch.jsds, "jsds/prochlorococcus.txt", sep = "\t")


# David -------------------------------------------------------------------

david <- fread(paste0(data.dir, "david/david.otus"),
               col.names = c("sample", "otu", "count")
               )
# parse subject and timepoints
david[, c("subject", "day") := tstrsplit(sample, "_")]
david[, rel.abund := count / sum(count), by = .(subject, day)]
david.jsds <- david[, jsds_from_d(.SD, "sample", "otu", "rel.abund")]
david.jsds[, c("subject.i", "day.i") := tstrsplit(sample.i, "_")]
david.jsds[, c("subject.j", "day.j") := tstrsplit(sample.j, "_")]
david.jsds[, c("sample.i", "sample.j") := NULL]

fwrite(david.jsds, "jsds/david.txt", sep = "\t")


# Gordon (cholera) --------------------------------------------------------


