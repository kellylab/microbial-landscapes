# setup -------------------------------------------------------------------


library(data.table)
library(philentropy)
library(tidyverse)

scripts.dir <- "./"
data.dir <- "../data/"
out.dir <- "jsds/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}


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
jsd <- as.data.table(jsd, keep.rownames = "sample.x")
jsd <- melt(jsd, id.vars = "sample.x", variable.name= "sample.y",
            value.name = "jsd")

fwrite(jsd, paste0(out.dir, "/cholera.txt"), sep = "\t")


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
jsd <- as.data.table(jsd, keep.rownames = "sample.x")
jsd <- melt(jsd, id.vars = "sample.x", variable.name= "sample.y",
            value.name = "jsd")
fwrite(jsd, paste0(out.dir, "/david.txt"), sep = "\t")


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
jsds <- as.data.table(jsds, keep.rownames = "sample.x")
jsds <- melt(jsds, id.vars = "sample.x", variable.name= "sample.y",
            value.name = "jsd")
fwrite(jsds, paste0(out.dir, "/prochlorococcus.txt"), sep = "\t")
