# setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(tidygraph)
library(cowplot)


jsd.dir <- "jsds/"
figs.dir <- "../../figures/tda/"
output.dir <- "../../output/mapper/"
if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
scripts.dir <- "../r/"
for (script in list.files("utils/", full.names = TRUE)) source(script)

write.graph <- function(tbl_graph, v2p, directory) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  tbl_graph %>%
    activate(nodes) %>%
    as.data.table %>%
    fwrite(paste0(directory, "vertices.txt"), sep = "\t", na = "NA")
  tbl_graph %>%
    activate(edges) %>%
    as.data.table %>%
    fwrite(paste0(directory, "edges.txt"), sep = "\t", na = "NA")
  names(v2p) <- c("point", "vertex")
  fwrite(v2p, paste0(directory, "/vertices-to-points.txt"), sep = "\t")
}

# cholera -----------------------------------------------------------------

source("../r/load_cholera_data.R")

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
ni <- c(15, 15)
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
theme_set(theme_graph(base_family = "Helvetica"))
ggraph(graf, "fr", niter = 1000) +
  geom_edge_link0() +
  geom_node_point(aes(color = f.state)) +
  scale_color_distiller(palette = "Spectral") +
  coord_equal()
graf <- graf %>%
  mutate(membership = components(.)$membership) %>%
  mutate(in.singleton = in.singleton(v2p$point.name, v2p$vertex, membership))
# assign minima and basins
graf <- graf %>% mutate(scaled.knn = mean.knn / size)
graf <- assign.basins(graf, "scaled.knn", ignore.singletons = TRUE)
graf <- graf %>%
  mutate(basin = mapply(function(b, x) if (x) NA else b,
                        b = basin, x = in.singleton)) %>%
  mutate(is.extremum = mapply(function(b, x) if (x) NA else b,
                        b = is.extremum, x = in.singleton))
write.graph(graf, v2p[, .(point.name, vertex)], paste0(output.dir, "cholera/"))


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
js.dist <- reshape2::acast(jsds, sample.x ~ sample.y, value.var = "jsd") %>%
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
ni <- c(30, 30)
po <- 50
mpr <- mapper2D(js.dist, ftr, num_intervals = ni, percent_overlap = po)
v2p <- vertex.2.points(mpr$points_in_vertex)
v2p <- merge(v2p, samples, by.x = "point.name", by.y = "sample")
vertices <- v2p[, .(
  f.subject = sum(subject == "A") / .N,
  mean.knn = mean(kNN),
  size = .N
), by = .(vertex, vertex.name)]

graf <- mapper.2.igraph(mpr)
graf <- graf %>%
  as_tbl_graph %>%
  activate(nodes) %>%
  left_join(vertices, by = c("name" = "vertex.name"))
graf <- graf %>%
  mutate(membership = components(.)$membership) %>%
  mutate(in.singleton = in.singleton(v2p$point.name, v2p$vertex, membership))

#' Find basins of attraction:
graf <- graf %>%
  activate(nodes) %>%
  mutate(scaled.knn = mean.knn / size)
graf <- assign.basins(graf, "scaled.knn", ignore.singletons = TRUE,
                      giant.only = FALSE) %>%
  activate(nodes) %>%
  mutate(basin = mapply(function(b, x) if (x) NA else b,
                        b = basin, x = in.singleton)) %>%
  mutate(is.extremum = mapply(function(b, x) if (x) NA else b,
                        b = is.extremum, x = in.singleton))
write.graph(graf, v2p[, .(point.name, vertex)], paste0(output.dir, "david/"))

# prochloroccoccus --------------------------------------------------------


month.2.phase <- function(month) {
  (month - 3) %% 12 / 12
}

source(paste0(scripts.dir, "load-prochlorococcus-data.R"))
jsds <- fread(paste0(jsd.dir, "prochlorococcus.txt"))
jsds[, distance := sqrt(jsd)]
dist.mat <- reshape2::acast(jsds, sample.x ~ sample.y, value.var = "distance")
samples <- unique(prochlorococcus[, -c("ecotype", "abundance")])
samples[, phase := month.2.phase(cal.month)]
k <- floor(nrow(samples) / 10)
knn <- dist2knn(dist.mat, k)
hist(knn)
set(samples, NULL, "knn", knn[samples$sample])
mds <- cmdscale(dist.mat, eig = TRUE) # 2D

#' GOF is good, but we find that samples are very unevenly distributed across
#' the 2D MDS-space, which is bad for Mapper:
mds$GOF # not too bad actually
plot(mds$points)

#' Converting to rank alleviates the problem somewhat.
#' Marginal distributions will be uniform, by definition.
rk.mds <- apply(mds$points, 2, rank, ties.method = "first")
plot(rk.mds)

#' Mapper call:
po <- 60
ni <- c(20, 20)
nb <- 10
ftr <- list(rk.mds[, 1], rk.mds[, 2])
mpr <- mapper2D(dist.mat, ftr,
                percent_overlap = po, num_intervals = ni,
                num_bins_when_clustering = nb)
v2p <- vertex.2.points(mpr$points_in_vertex)
v2p$sample <- rownames(dist.mat)[v2p$point]
setkey(v2p, sample)
setkey(samples, sample)
v2p <- samples[v2p]
v2p[, depth := as.numeric(depth)]
vertices <- v2p[, .(size = .N,
                    f.bats = sum(site == "bats") / .N,
                    mean.depth = mean(depth, na.rm = TRUE),
                    mean.temp = mean(temp, na.rm = TRUE),
                    mean.sal = mean(sal, na.rm = TRUE),
                    mean.calmonth = mean(cal.month, na.rm = TRUE),
                    median.depth = median(depth, na.rm = TRUE),
                    median.temp = median(temp, na.rm = TRUE),
                    median.sal = median(sal, na.rm = TRUE),
                    mean.knn = mean(knn, na.rm = TRUE)
                    ),
                by = .(vertex, vertex.name)]
graf <- mapper.2.igraph(mpr) %>%
  as_tbl_graph %>%
  left_join(vertices, by = c("name" = "vertex.name")) %>%
  mutate(scaled.knn = mean.knn / size)
graf <- assign.basins(graf, "scaled.knn", ignore.singletons = TRUE)
graf <- graf %>%
  mutate(membership = components(.)$membership) %>%
  mutate(in.singleton = in.singleton(v2p$point.name, v2p$vertex, membership))
graf <- mutate(graf, basin = nullify(in.singleton, basin))
write.graph(graf, v2p[, .(point.name, vertex)],
            paste0(output.dir, "prochlorococcus/"))


# nahant ------------------------------------------------------------------


source(paste0(scripts.dir, "load-nahant-data.R"))

nahant <- nahant[!is.na(kingdom)]
nahant[, freq := value / sum(value), by = .(kingdom, day)]

# jsds
jsds <- c(Bacteria = "Bacteria", Eukaryota = "Eukaryota") %>%
  lapply(function(k) {
    fn <- paste0("jsds/nahant-", k, ".txt")
    d <- fread(fn)
    m <- dcast(d, day.i ~ day.j, value.var = "jsd")
    rn <- m$day.i
    m <- as.matrix(m[, -1])
    rownames(m) <- rn
    m
})
distances <- lapply(jsds, sqrt)

#' 2D MDS is lossy, but suggests existence of 2 bacterial and 2-3 eukaryotic
#' states:
mds <- lapply(distances, cmdscale, eig = TRUE)
for (x in mds) {
  print(x$GOF)
  plot(x$points, asp = 1)
}
rk.mds <- lapply(mds, function(x) {
  pts <- x$points
  apply(pts, 2, rank, ties.method = "first")
})
for (x in rk.mds) {
  plot(x, asp = 1)
}

# mapper call
bn <- 10
en <- 5
ni <- list(c(bn, bn), c(en, en))
po <- list(70, 70)
bb <- 10
eb <- 20
nbin <- list(bb, eb)
mpr <- mapply(function(distance, rk.mds, po, ni, nb) {
  mapper2D(distance, list(rk.mds[, 1], rk.mds[, 2]), num_intervals = ni,
           percent_overlap = po, num_bins_when_clustering = nb)
}, distance = distances, rk.mds = rk.mds, ni = ni, po = po, nb = nbin,
SIMPLIFY = FALSE)
summaries <- lapply(mpr, summary, plot = TRUE)
plot_grid(plotlist = summaries, ncol = 2, labels = names(summaries), align = "h")

# make graphs
grafs <- lapply(mpr, mapper.2.igraph) %>%
  lapply(as_tbl_graph)

# map vertices to days
v2p <- lapply(mpr, function(x) {
  dt <- vertex.2.points(x$points_in_vertex)
  setnames(dt, "point.name", "day")
  dt[, day := as.numeric(day)]
  dt
})

# knn
knn <- lapply(distances, function(d) {
  k <- round(nrow(d) / 10)
  v <- apply(d, 1, k.first, k = k)
  data.table(day = as.numeric(names(v)), kNN = v)
})
v2p <- mapply(merge, x = v2p, y = knn, MoreArgs = list(by = "day"),
              SIMPLIFY = FALSE)
vertices <- lapply(v2p, function(dt) {
  dt[, .(mean.day = mean(day), mean.knn = mean(kNN), size = .N), by = vertex]
})
grafs <- mapply(function(g, d) {
  d$name <- paste0("v", d$vertex)
  g %>%
    activate(nodes) %>%
    left_join(d, by = "name") %>%
    mutate(scaled.knn = mean.knn / size)
}, g = grafs, d = vertices, SIMPLIFY = FALSE)
grafs <- lapply(grafs, assign.basins, fn = "scaled.knn", ignore.singletons = TRUE)
grafs <- lapply(grafs, function(g) {
  g %>% activate(nodes) %>% mutate(basin = as.factor(basin))
})

# basin plot
set.seed(0)
pl <- lapply(grafs, function(g) {
  ggraph(g, "fr") +
    geom_edge_link0() +
    geom_node_point(aes(size = size, fill = basin), shape = 21) +
    theme_graph(base_family = "Helvetica") +
    theme(aspect.ratio = 1)
})
plot_grid(plotlist = pl, labels = names(pl))

# save
dirs <- list(paste0(output.dir, "nahant/bacteria/"),
             paste0(output.dir, "nahant/eukaryotes/"))
mapply(write.graph, tbl_graph = grafs,
       v2p = lapply(v2p, function(d) d[, .(point.name = day, vertex)]),
       directory = dirs)
