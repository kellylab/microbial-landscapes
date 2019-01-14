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

# validation parameters
nrep <- 10             # number of replicates
rs <- c(0.9, 0.5, 0.1) # downsampling ratios

jsd.dir <- "jsds/"
figs.dir <- "../../figures/tda/"
output.dir <- "../../output/mapper/"
if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
scripts.dir <- "../r/"
for (script in list.files("utils/", full.names = TRUE)) source(script)

#' Construct function to call Mapper2D and map vertices to original data points
#'
#' @param ni      number intervals
#' @param po      percent overlap
#' @param vfn     optional, function for aggregating per-vertex data
#'
#' @return list, [1] Mapper graph [2] vertex-point map
#' @export
#'
#' @examples
mapper2.call <- function(ni, po, vfn = NULL) {
  function(dist, samples, ftr) {
    mpr <- mapper2D(dist, filter_values = ftr, num_intervals = ni,
                    percent_overlap = po)
    v2p <- vertex.2.points(mpr$points_in_vertex)
    v2p <- merge(v2p, samples, by.x = "point.name", by.y = "sample")
    v2p[, frac := 1 / .N, by = point]
    graf <- mapper.2.igraph(mpr) %>%
      as_tbl_graph
    graf <- graf %>%
      mutate(membership = components(.)$membership) %>%
      mutate(in.singleton = in.singleton(v2p$point.name, v2p$vertex, membership))
    if (!is.null(vfn)) {
      vertices <- do.call(vfn, list(map = v2p))
      setorder(vertices, vertex)
      # browser()
      graf <- graf %>%
        activate(nodes) %>%
        left_join(vertices, by = c("name" = "vertex.name"))
    }
    list(graph = graf, map = v2p)
  }
}

#' Random subset of distance matrix and data
#'
#' @param dist distance matrix
#' @param dt   original data
#' @param r    fraction to sample
#'
#' @return list, [1] reduced distance matrix [2] reduced data
#' @export
#'
#' @examples
subsample <- function(dist, dt, r = 0.9) {
  n <- nrow(dist)
  size <- round(n * r)
  idx <- sample.int(n, size)
  dist <- dist[idx, idx]
  samps <- rownames(dist)
  setkey(dt, sample)
  dt <- dt[samps]
  list(dist = dist, data = dt)
}

#' Validate Mapper representation by subsampling by given ratios some number of
#' times
#'
#' @param dist  original distance matrix
#' @param dt    original sample data
#' @param fn    function to wrap Mapper
#' @param rs    vector of downsampling coefficients
#' @param nrep  number of replicates per coefficient
#'
#' @return list of list of Mapper representations and vertex-sample maps
#' @export
#'
#' @examples
validate <- function(dist, dt, fn, rs, nrep) {
  lapply(rs, function(r, nrep, dist, dt, fn) {
    lapply(seq_len(nrep), function(i, r, dist, dt, fn) {
      subsamp <- subsample(dist, dt, r)
      mds <- cmdscale(subsamp$dist)
      rk.mds <- apply(mds, 2, rank)
      do.call(fn, list(dist = subsamp$dist,
                       samples = subsamp$data,
                       ftr = list(rk.mds[,1], rk.mds[,2])))
    }, r, dist, dt, fn)
  }, nrep, dist, dt, fn)
}

#' Basic FR layout and plot of Mapper graph
#'
#' @param graf
#' @param edge
#' @param node
#' @param exclude.singletons
#' @param seed
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot.mapper.graph <- function(graf,
                              edge = geom_edge_link0(),
                              node = geom_node_point(),
                              exclude.singletons = FALSE,
                              seed = NULL) {
  theme_set(theme_graph(base_family = "Helvetica"))
  set.seed(seed)
  if (exclude.singletons) {
    graf <- graf %>%
      activate(nodes) %>%
      filter(!in.singleton)
  }
  ggraph(graf, "fr", niter = 1000) +
    edge +
    node +
    theme_graph(base_family = "Helvetica") +
    theme(aspect.ratio = 1)
}

#' Plot downsampled Mapper graphs in a grid
#'
#' @param subsets  list of lists of downsampled Mappers
#' @param fn       function used for plotting each downsampled Mapper
#'
#' @return
#' @export         list of ggplot objects
#'
#' @examples
batch.plot <- function(subsets, fn) {
  pl <- lapply(subsets, function(l) {
    plots <- lapply(l, function(mpr) {
      do.call(fn, list(graf = mpr$graph)) +
        theme(legend.position = "none",
              plot.margin = unit(c(0, 0, 0, 0), "points"))
    })
    plot_grid(plotlist = plots, nrow = 2)
  })
}

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
ftr <- list(rk.mds[,1], rk.mds[,2])

cholera.vertices <- function(map) {
  map[, .(mean.knn = mean(knn),
          f.state = sum(diagnosis == "diarrhea") / .N,
          mean.t = mean(hour),
          size = .N),
      by = .(vertex, vertex.name)]
}

cholera.mapper <- mapper2.call(ni, po, cholera.vertices)
mpr <- cholera.mapper(js.dist, gordon.samples, ftr)

plot.fstate <- function(graf, ...) {
  plot.mapper.graph(graf,
                    node = geom_node_point(aes(color = f.state, size = size)),
                    ...) +
    scale_color_distiller(palette = "Spectral")
}

plot.fstate(mpr$graph, seed = 1)

# validate
subsets <- validate(js.dist, gordon.samples, cholera.mapper, rs, nrep)
pl <- batch.plot(subsets,
                 function(graf) plot.fstate(graf) + scale_size_area(max_size = 2))
plot_grid(plotlist = pl, ncol = 1, labels = rs)

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

david.vertices <- function(map) {
  map[, .(f.subject = sum(subject == "A") / .N,
          mean.knn = mean(kNN),
          size = .N),
      by = .(vertex, vertex.name)]
}
david.mapper <- mapper2.call(ni, po, david.vertices)
mpr <- david.mapper(js.dist, david.samples, ftr)

plot.fsubject <- function(graf) {
  plot.mapper.graph(graf,
                    node = geom_node_point(aes(color = f.subject, size = size))) +
    scale_color_distiller(palette = "Spectral")
}

# validate
subsets <- validate(js.dist, david.samples, david.mapper, rs, nrep)
pl <- batch.plot(subsets, function(graf) {
  plot.fsubject(graf) + scale_size_area(max_size = 2)
})

plot_grid(plotlist = pl, ncol = 1, labels = rs)

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
proc.vertices <- function(map) {
  map[, .(size = .N,
          f.bats = sum(site == "bats") / .N,
          mean.depth = mean(depth, na.rm = TRUE),
          mean.temp = mean(temp, na.rm = TRUE),
          mean.sal = mean(sal, na.rm = TRUE),
          mean.calmonth = mean(cal.month, na.rm = TRUE)
          ),
      by = .(vertex, vertex.name)]
}
proc.mapper <- mapper2.call(ni, po, proc.vertices)
# nb <- 10
ftr <- list(rk.mds[, 1], rk.mds[, 2])
mpr <- proc.mapper(dist.mat, prochlorococcus.samples, ftr)

plot.depth <- function(graf) {
  plot.mapper.graph(graf,
                    node = geom_node_point(aes(color = mean.depth, size = size)),
                    exclude.singletons = TRUE,
                    seed = 0) +
    scale_color_distiller(palette = "Blues", direction = 1)
}

plot.depth(mpr$graph)

# validate
subsets <- validate(dist.mat, prochlorococcus.samples, proc.mapper, rs, nrep)
pl <- batch.plot(subsets, function(graf) {
  plot.depth(graf) + scale_size_area(max_size = 2)
})
plot_grid(plotlist = pl, ncol = 1, labels = rs)

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

# summary statistics for vertices
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
set.seed(0)
pl <- lapply(grafs, function(g) {
  ggraph(g, "fr") +
    geom_edge_link0() +
    geom_node_point(aes(size = size, fill = mean.day), shape = 21) +
    scale_fill_distiller(palette = "Spectral") +
    theme_graph(base_family = "Helvetica") +
    theme(aspect.ratio = 1)
})
plot_grid(plotlist = pl, labels = names(pl))

# knn and basins
grafs <- lapply(grafs, assign.basins, fn = "scaled.knn", ignore.singletons = TRUE)
grafs <- lapply(grafs, function(g) {
  g %>% activate(nodes) %>% mutate(basin = as.factor(basin))
})
grafs <- mapply(function(graf, v2p) {
  comps <- components(graf)
  graf %>%
    activate(nodes) %>%
    mutate(in.singleton = in.singleton(v2p$day, v2p$vertex, comps$membership)) %>%
    mutate(basin = nullify(in.singleton, basin)) %>%
    mutate(basin = as.factor(basin))
}, graf = grafs, v2p = v2p, SIMPLIFY = FALSE)
set.seed(0)
layouts <- lapply(grafs, create_layout, layout = "fr")
pl <- mapply(function(g, lo) {
  ggraph(g, "manual", node.positions = lo[, c("x", "y")]) +
    geom_edge_link0() +
    geom_node_point(aes(size = size, fill = basin), shape = 21) +
    theme_graph(base_family = "Helvetica") +
    theme(aspect.ratio = 1)
}, g = grafs, lo = layouts, SIMPLIFY = FALSE)
plot_grid(plotlist = pl, labels = names(pl))

# save
dirs <- list(paste0(output.dir, "nahant/bacteria/"),
             paste0(output.dir, "nahant/eukaryotes/"))
mapply(write.graph, tbl_graph = grafs,
       v2p = lapply(v2p, function(d) d[, .(point.name = day, vertex)]),
       directory = dirs)

# sketch for downstream analyses
for (d in v2p) setkey(d, day)
map <- merge(v2p$Bacteria, v2p$Eukaryota, suffixes = paste0(".", names(v2p)),
             allow.cartesian = TRUE) %>%
  .[, .(vertex.Bacteria = paste0("Bacteria", vertex.Bacteria),
        vertex.Eukaryota = paste0("Eukaryota", vertex.Eukaryota))] %>%
  unique

# rename vertices to merge graphs
grafs <- mapply(function(g, k) {
  mutate(g, name = paste0(k, vertex), basin = paste0(k, basin), type = k)
}, g = grafs, k = names(grafs), SIMPLIFY = FALSE)

# offset the euk vertices
delx <- 2 * max(layouts$Bacteria$x) - min(layouts$Bacteria$x)
merged.layout <- layouts$Eukaryota %>%
  mutate(x = x + delx) %>%
  rbind(layouts$Bacteria)
merged.graf <- graph_join(grafs$Eukaryota, grafs$Bacteria)
merged.graf <- merged.graf %>%
  activate(edges) %>%
  mutate(original = TRUE)
for (i in seq(nrow(map))) {
  merged.graf <- merged.graf +
    edge(map[i,]$vertex.Bacteria, map[i,]$vertex.Eukaryota, original = FALSE)
}
merged.graf <- as_tbl_graph(merged.graf)
merged.graf <- activate(merged.graf, nodes) %>%
  arrange(type, basin, -size)
ggraph(merged.graf, "linear", circular = TRUE) +
  geom_edge_link0(data = function(d) filter(get_edges()(d), original == FALSE),
                  alpha = 0.1) +
  geom_node_point(aes(size = size, shape = type, color = basin)) +
  theme_graph(base_family = "Helvetica") +
  theme(aspect.ratio = 1)

# connectivity table btwn bac and euk basins
graf.df <- lapply(grafs, as.data.frame)
map <- merge(map, graf.df$Bacteria, by.x = "vertex.Bacteria", by.y = "name")
map <- merge(map, graf.df$Eukaryota, by.x = "vertex.Eukaryota", by.y = "name",
             suffixes = paste0(".", names(grafs)))
ggplot(map, aes(x = basin.Bacteria, y = basin.Eukaryota)) +
  stat_bin_2d() +
  coord_equal()
