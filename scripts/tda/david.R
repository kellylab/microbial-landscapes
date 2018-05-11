library(philentropy)
library(TDAmapper)
library(data.table)
library(tidyverse)
library(ggraph)
library(igraph)
library(cowplot)

# setup -------------------------------------------------------------------


scripts.dir <- "../r/"
source(paste0(scripts.dir, "load_david_data.R"))
samples <- unique(david[, .(sample, subject, day)])
samples[, day := as.numeric(day)]

# if haven't pre-computed jsds, run make-jsds-david.R
jsds <- fread("jsds/david.txt")
jsds[, distance := sqrt(jsd)]
dist.mat <- dcast(jsds, sample.x ~ sample.y, value.var = "jsd")
sample.names <- dist.mat$sample.x
dist.mat <- as.matrix(dist.mat[, -1])
rownames(dist.mat) <- sample.names

#' # Preview of data using 2D MDS
#'
#' 2D MDS shows much more variance in dimension 1 than dimension 2, resulting in
#' large empty spaces.
mds2d <- as.data.table(cmdscale(dist.mat, 2), keep.rownames = "sample")
setnames(mds2d, c("V1", "V2"), c("mds1", "mds2"))
mds2d[, subject := tstrsplit(sample, "_")[[1]]]
ggplot(mds2d, aes(x = mds1, y = mds2)) +
  geom_point(aes(color = subject))

#' # Mapper
#' ## Filter by 2D MDS
# mds filter -------------------------------------------------------------
mpr <- mapper1D(dist.mat, cmdscale(dist.mat, 1), num_intervals = 100,
                percent_overlap = 80)
mpr$points_in_vertex <- lapply(mpr$points_in_vertex, function(v, rn) {
  names(v) <- rn[v]
  v
}, rn = rownames(dist.mat))

# mapper output ----------------------------------------------------


source("mapper.adj.R")
mpr$adjacency <- mapper.adj(mpr$points_in_vertex) # correct adj matrix

#' Format Mapper output:
source("vertex.2.points.R")
v2p <- vertex.2.points(mpr$points_in_vertex) # map vertices to samples
v2p <- merge(v2p, samples, by.x = "point.name", by.y = "sample")
source("mapper.2.igraph.R")
graf <- mapper.2.igraph(mpr)
V(graf)$subject <- sapply(mpr$points_in_vertex, function(pts, dt) {
  pnames <- names(pts)
  setkey(dt, sample)
  dt[pnames, sum(subject == "A") / .N]
}, dt = samples)


#' ## Plots

source("plot.mapper.R")
#' Fraction of samples in each vertex belonging to each subject:
set.seed(0)
lo <- create_layout(graf, "fr")
# store xy positions in graph for later
V(graf)$x <- lo$x
V(graf)$y <- lo$y
plot.mapper(lo, aes_(size = ~size, color = ~subject),
            list(color = "fraction A")) +
  scale_color_gradient2(midpoint = 0.5, mid = "yellow")

#' ### Eric's trip
#'
# eric’s frames -----------------------------------------------------------


source("sample.subgraphs.R")
sample.sgrafs <- sample.subgraphs(
  v2p[, .(samples = point.name, vertices = vertex)],
  graf)
# eventmetadata
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
setkey(v2p, subject, day)
setkey(events, subject, day)
v2p <- events[v2p]
b.samples <- samples["B", sample]
setkey(samples, sample)
# setkey(v2p, point.name)
b.events <- events[samples[b.samples, .(subject, day)], event]
sample.overlays <- mapply(function(s, sg, e, lo) {
  if (!is.null(sg)) {
    if (grepl("pre", e)) {
      clr <- "blue"
    } else if (grepl("post", e)) {
      clr <- "yellow"
    } else {
      clr <- "red"
    }
    ttl <- paste(s, e, sep = "; ")
    g <- plot.mapper(lo, aes_(size = ~size), list(title = ttl),
                     color = "grey50")
    g + geom_point(aes(x = x, y = y, size = size),
                   data = as.data.frame(vertex_attr(sg)),
                   color = clr)



  } else {
    NULL
  }
}, s = b.samples, sg = sample.sgrafs[b.samples], e = b.events,
MoreArgs = list(lo = lo), SIMPLIFY = FALSE)
setorder(samples, day)
setkey(samples, subject)
sample.overlays <- sample.overlays[samples["B", sample]]
sample.overlays <- sample.overlays[sapply(sample.overlays,
                                          function(x) !is.null(x))]
for (i in seq_along(sample.overlays)) {
  so <- sample.overlays[[i]]
  save_plot(paste0("david-subject-trajectories/frames/B", i, ".png"), so)
}
#'
#'
#' source("vertex.attribute.R")
#'
#' subj.t <- function(vn, subj, v2p) {
#'   setkey(v2p, vertex.name, subject)
#'   v2p[.(vn, subj), mean(day)]
#' }
#' V(graf)$mean.tB <- sapply(V(graf)$name, subj.t, subj = "B", v2p = v2p)
#' set.seed(0)
#' lo <- create_layout(graf, "fr")
#' plot.mapper(lo, aes_(color = ~mean.tB), list(color = "B's mean day")) +
#'   scale_color_distiller(palette = "Spectral")
#'
#' #' ### Events
