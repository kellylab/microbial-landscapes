# setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(cowplot)


util.dir <- "../r/"
jsd.dir <- "jsds/"
figs.dir <- "../../figures/tda/"
output.dir <- "../../output/mapper/"

source("vertex.2.points.R")
source("dist2knn.R")
source("sample.subgraphs.R")
source("assign.basins.R")
source("mapper.2.igraph.R")

write.graph <- function(tbl_graph, directory) {
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
  tbl_graph %>%
    activate(nodes) %>%
    as.data.table %>%
    fwrite(paste0(output.dir, paste0(directory, "vertices.txt")),
                  sep = "\t", na = "NA")
  tbl_graph %>%
    activate(edges) %>%
    as.data.table %>%
    fwrite(paste0(output.dir, paste0(directory, "edges.txt")),
                  sep = "\t", na = "NA")
}

# cholera -----------------------------------------------------------------

source(paste0(util.dir, "load_cholera_data.R"))

gordon[, freq := count / sum(count), by = sample]
gordon.samples <- unique(gordon[, .(sample, subject, diagnosis, id, hour)])
gordon.samples[, idx := frank(hour), by = subject]
gordon.samples[, progression := idx / max(idx), by = subject]
distribs <- dcast(gordon, sample ~ otu, value.var = "freq", fill = 0)
sample.names <- distribs$sample
distribs <- as.matrix(distribs[, -1])

# import js distance
jsd.file <- paste0(jsd.dir, "cholera.txt")
if (!file.exists(jsd.file)) {
  print("Making cholera distance matrix...")
  source("make-jsds-cholera.R")
}
jsd <- fread(jsd.file)
js.dist <- dcast(jsd, sample.x ~ sample.y, value.var = "jsd") %>%
  column_to_rownames("sample.x") %>%
  data.matrix(rownames.force = TRUE) %>%
  sqrt

#' ## K-nearest neighbor density
k <- round(nrow(gordon.samples) / 10)
knn <- dist2knn(js.dist, k = k)
gordon.samples[, knn := knn[sample]]
ggplot(gordon.samples, aes(x = progression, y = knn)) +
  geom_smooth() +
  geom_point(aes(color = diagnosis)) +
  facet_wrap(~ subject)

#' Diarrhea samples only showing disease duration:
setkey(gordon.samples, diagnosis)
theme_set(theme_cowplot())
ggplot(gordon.samples["diarrhea"], aes(x = hour, y = knn)) +
  geom_smooth() +
  geom_point() +
  facet_wrap(~ subject)

#' MDS sketch suggests existence of at least 2 clusters.
mds2 <- cmdscale(js.dist, 2, eig = TRUE)
mds2$GOF # lossy
rk.mds <- apply(mds2$points, 2, rank)
plot(rk.mds)

# mapper call
po <- 70
ni <- c(10, 10)
mpr <- mapper2D(js.dist, filter_values = list(rk.mds[,1], rk.mds[,2]),
                percent_overlap = po,
                num_intervals = ni)
v2p <- vertex.2.points(mpr$points_in_vertex)
v2p <- merge(v2p, gordon.samples, by.x = "point.name", by.y = "sample")
v2p[, frac := 1 / .N, by = point]
vertices <- v2p[, .(mean.knn = mean(knn),
                    f.state = sum(diagnosis == "diarrhea") / .N,
                    mean.t = mean(hour),
                    size = .N),
                by = .(vertex, vertex.name)]
setorder(vertices, vertex)
graf <- mapper.2.igraph(mpr) %>%
  as_tbl_graph %>%
  activate(nodes) %>%
  left_join(vertices, by = c("name" = "vertex.name"))
# assign minima and basins
graf <- assign.basins(graf, "mean.knn", ignore.singletons = TRUE)
write.graph(graf, paste0(output.dir, "cholera/"))


# david -------------------------------------------------------------------


source(paste0(scripts.dir, "load_david_data.R"))
samples <- unique(david[, .(sample, subject, day)])
samples[, day := as.numeric(day)]
jsd.file <- paste0(jsd.dir, "david.txt")
if (!file.exists(jsd.file)) {
  print("Making David sample JSDs...")
  source("make-jsds-david.R")
}
jsds <- fread(jsd.file)
js.dist <- dcast(jsds, sample.x ~ sample.y, value.var = "jsd") %>%
  column_to_rownames("sample.x") %>%
  data.matrix(rownames.force = TRUE) %>%
  sqrt
k <- round(nrow(samples) / 10) # calculate k nearest neighbors
kNN <- dist2knn(js.dist, k)
samples[, kNN := kNN[sample]]

# event metadata
setkey(samples, subject)
events <- mapply(function(start, end, subject, event) {
  x <- data.table(day = seq(start, end, by = 1), subject = subject,
                  event = event)
  x
}, start = c(0, 71, 80, 104, 123, 0, 151, 160 ),
end = c(70, 122, 85, 113, samples["A", max(day)],
          150, 159, samples["B", max(day)]),
subject = c("A", "A", "A", "A", "A",
            "B", "B", "B"),
event = c("US (pre)", "travel", "diarrhea 1", "diarrhea 2", "US (post)",
          "pre-Salmonella", "Salmonella", "post-Salmonella"),
SIMPLIFY = FALSE) %>% rbindlist(use.names = TRUE)
# collapse event labels per day
events <- events[, .(event = paste(event, collapse = " + ")),
                 by = .(subject, day)]
samples <- merge(samples, events, by = c("subject", "day"))

#' 2D MDS shows uneven distribution. Transform to rank to make more even.
mds2d <- cmdscale(js.dist, eig = TRUE)
mds2d$GOF # not very good
plot(mds2d$points)
rk.mds <- apply(mds2d$points, 2, rank, ties.method = "first")
plot(rk.mds)

# mapper call

ftr <- list(rk.mds[, 1], rk.mds[, 2])
ni <- c(20, 20)
po <- 70
mpr <- mapper2D(js.dist, ftr, num_intervals = ni, percent_overlap = po)
v2p <- vertex.2.points(mpr$points_in_vertex)
v2p <- merge(v2p, samples, by.x = "point.name", by.y = "sample")
vertices <- v2p[, .(
  subject = sum(subject == "A") / .N,
  mean.knn = mean(kNN),
  size = .N
), by = .(vertex, vertex.name)]

graf <- mapper.2.igraph(mpr)
graf <- graf %>%
  as_tbl_graph %>%
  activate(nodes) %>%
  left_join(vertices, by = c("name" = "vertex.name"))

#' Find basins of attraction:
graf <- assign.basins(graf, "mean.knn", ignore.singletons = TRUE)
write.graph(graf, paste0(output.dir, "david/"))

# prochloroccoccus --------------------------------------------------------


