# setup -------------------------------------------------------------------


library(data.table)
library(philentropy)
library(tidyverse)

scripts.dir <- "../r/"


# cholera ---------------------------------------------------------------

source(paste0(scripts.dir, "load_cholera_data.R"))

gordon[, freq := count / sum(count), by = sample]
gordon.samples <- unique(gordon[, .(sample, subject, diagnosis, id, hour)])
gordon.samples[, idx := frank(hour), by = subject]
distribs <- dcast(gordon, sample ~ otu, value.var = "freq", fill = 0)
sample.names <- distribs$sample
distribs <- as.matrix(distribs[, -1])

# compute js distance
jsd <- JSD(distribs)
rownames(jsd) <- sample.names
colnames(jsd) <- sample.names
jsd <- melt(jsd, varnames = c("sample.x", "sample.y"),
            value.name = "jsd")

fwrite(jsd, "jsds/cholera.txt", sep = "\t")


# david -------------------------------------------------------------------


source(paste0(scripts.dir, "load_david_data.R"))
david[, rel.abundance := count / sum(count), by = sample]
sample.abundances <- dcast(david, sample ~ otu, value.var = "rel.abundance",
                           fill = 0)
sample.names <- sample.abundances$sample
sample.abundances <- as.matrix(sample.abundances[, -1])

# pairwise distances
jsd <- JSD(sample.abundances)
rownames(jsd) <- sample.names
colnames(jsd) <- sample.names

# format and write
jsd <- melt(jsd, varnames = c("sample.x", "sample.y"), value.name = "jsd")
fwrite(jsd, "jsds/david.txt", sep = "\t")


# prochlorococcus ---------------------------------------------------------


source(paste0(scripts.dir, "load-prochlorococcus-data.R"))

prochlorococcus[, sample := paste(site, cruiseid, depth, sep = "-")]
prochlorococcus[, total.abundance := sum(abundance), by = sample]
prochlorococcus <- prochlorococcus[total.abundance > 0]
prochlorococcus[, freq := abundance / total.abundance, by = sample]
distribs <- dcast(prochlorococcus, sample ~ ecotype, value.var = "freq")
samples <- distribs[, sample]
jsds <- JSD(as.matrix(distribs[, -1]))
rownames(jsds) <- samples
colnames(jsds) <- samples
jsds <- as.data.table(melt(jsds, varnames = c("sample.x", "sample.y"),
                           value.name = "jsd"))
fwrite(jsds, "jsds/prochlorococcus.txt", sep = "\t")


# nahant ------------------------------------------------------------------


source(paste0(scripts.dir, "load-nahant-data.R"))

nahant <- nahant[!is.na(kingdom)]
nahant[, freq := value / sum(value), by = .(kingdom, day)]

# jsds
jsds <- lapply(split(nahant, by = "kingdom"), function(dt) {
  x <- dcast(dt, day ~ OTU, value.var = "freq", fill = 0)
  rn <- x$day
  x <- as.matrix(x[, -1])
  rownames(x) <- rn
  m <- JSD(x)
  rownames(m) <- rn
  colnames(m) <- rn
  m
})

jsds <- lapply(jsds, function(m) {
  d <- as.data.table(m, keep.rownames = "day.i")
  d <- melt(d, id.vars = "day.i", variable.name = "day.j", value.name = "jsd")
  d
})

for (kingdom in names(jsds)) {
  fn <- paste0("jsds/nahant-", kingdom, ".txt")
  fwrite(jsds[[kingdom]], fn, sep = "\t")
}