library(tidyverse)
library(data.table)
library(philentropy)
library(infograf)

source("../../../../wkc.r-utils/df_2_mat.R")
source("../../../../wkc.r-utils/mat_2_df.R")

# general function to get JSDs from pairwise data
jsds_from_d <- function(d, sample.cn, var.cn, val.cn, suf = c(".i", ".j")) {
  form <- paste(sample.cn, "~", var.cn)
  d <- dcast(d, form, value.var = val.cn)
  rn <- d[[sample.cn]]
  d <- select(d, -matches(sample.cn))
  m <- as.matrix(d)
  rownames(m) <- rn
  jsds <- JSD(m) # WILL BE NaN IF ANY ECOTYPE IS 0
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
# philentropy::JSD gives NaN if there are zeros
# find these pairs and calculate the JSD with my own function
setkey(prochlorococcus, site, sample)
proch.jsds[is.nan(jsd), jsd := {
  sit <- site
  x <- prochlorococcus[.(sit, c(sample.i, sample.j))]
  x <- dcast(x, sample ~ ecotype, value.var = "rel.abund")
	x <- as.data.frame(x)
  m <- df_2_mat(x)
  genJSD(t(m))
},
by = .(site, sample.i, sample.j)]
proch.jsds[, c("year.i", "month.i", "depth.i") := tstrsplit(sample.i, "_")]
proch.jsds[, c("year.j", "month.j", "depth.j") := tstrsplit(sample.j, "_")]
proch.jsds[, c("sample.i", "sample.j") := NULL]
fwrite(proch.jsds, "jsds/prochlorococcus.txt", sep = "\t")
