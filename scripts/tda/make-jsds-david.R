library(philentropy)
library(data.table)

# setup -------------------------------------------------------------------


scripts.path <- "../r/"
source(paste0(scripts.path, "load_david_data.R"))
david[, rel.abundance := count / sum(count), by = sample]
sample.abundances <- dcast(david, sample ~ otu, value.var = "rel.abundance",
                           fill = 0)
sample.names <- sample.abundances$sample
sample.abundances <- as.matrix(sample.abundances[, -1])

# pairwise distances ------------------------------------------------------


jsd <- JSD(sample.abundances)
rownames(jsd) <- sample.names
colnames(jsd) <- sample.names

# format and write
jsd <- melt(jsd, varnames = c("sample.x", "sample.y"), value.name = "jsd")
fwrite(jsd, "jsds/david.txt", sep = "\t")