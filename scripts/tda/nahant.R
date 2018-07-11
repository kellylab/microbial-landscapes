library(tidyverse)
library(data.table)
library(philentropy)
library(TDAmapper)
library(ggraph)
library(igraph)
library(tidygraph)
library(cowplot)

scripts.dir <- "../r/"
figs.dir <- "../../figures/tda/"
source(paste0(scripts.dir, "load-nahant-data.R"))
source("k.first.R")
source("vertex.2.points.R")
source("mapper.2.igraph.R")
source("plot.mapper.R")

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
k <- 10
samples <- lapply(distances, function(m, k) {
  v <- apply(m, 1, k.first, k = k)
  data.table(day = as.numeric(names(v)), kNN = v)
}, k = k)
#' Bacteria show trend of increasing instability toward day 250, while
#' eukaryotes have about the same stability:
ggplot(rbindlist(samples, idcol = "kingdom"), aes(x = day, y = kNN)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ kingdom, nrow = 2)
mds <- lapply(distances, cmdscale, eig = TRUE)
#' 2D MDS is lossy, but suggests existence of 2 bacterial and 2-3 eukaryotic
#' states:
for (x in mds) {
  print(x$GOF)
  plot(x$points)
}
rk.mds <- lapply(mds, function(x) {
  pts <- x$points
  apply(pts, 2, rank, ties.method = "first")
})


# mapper ------------------------------------------------------------------

ni <- c(10, 10)
po <- 70
mpr <- mapply(function(distance, rk.mds) {
  mapper2D(distance, list(rk.mds[, 1], rk.mds[, 2]), num_intervals = ni,
           percent_overlap = po)
}, distance = distances, rk.mds = rk.mds, SIMPLIFY = FALSE)
v2p <- lapply(mpr, function(x) {
  dt <- vertex.2.points(x$points_in_vertex)
  setnames(dt, "point.name", "day")
  dt[, day := as.numeric(day)]
  dt
})
v2p <- mapply(merge, x = v2p, y = samples, MoreArgs = list(by = "day"),
              SIMPLIFY = FALSE)
vertices <- lapply(v2p, function(dt) {
  dt[, .(mean.day = mean(day),
         mean.kNN = mean(kNN),
         size = .N),
     by = .(vertex, vertex.name)]
})
grafs <- lapply(mpr, mapper.2.igraph)
grafs <- mapply(function(graf, verts) {
  graf %>%
    as_tbl_graph %>%
    activate(nodes) %>%
    left_join(verts, by = c("name" = "vertex.name"))
}, graf = grafs, verts = vertices, SIMPLIFY = FALSE)
set.seed(0)
layouts <- lapply(grafs, create_layout, layout = "fr")
grafs <- mapply(function(g, l) {
  g %>%
    activate(nodes) %>%
    mutate(x = l$x, y = l$y)
}, g = grafs, l = layouts, SIMPLIFY = FALSE)
size.plots <- mapply(function(lo, kingdom) {
  p <- plot.mapper(lo, aes_(size = ~size, color = ~mean.day))
  p + scale_color_distiller(palette = "Spectral")
}, lo = layouts, kingdom = names(layouts), SIMPLIFY = FALSE)
plot_grid(plotlist = size.plots, labels = names(size.plots))
save_plot(paste0(figs.dir, "nahant-mean-day.pdf"), last_plot(),
          ncol = 2, base_height = 6)
kNN.plots <- mapply(function(lo, kingdom) {
  p <- plot.mapper(lo, aes_(size = ~size, color = ~mean.kNN),
                   list(title = kingdom))
  p + scale_color_distiller(palette = "Spectral") +
    theme(legend.position = "bottom")
}, lo = layouts, kingdom = names(layouts), SIMPLIFY = FALSE)
plot_grid(plotlist = kNN.plots)

#' # Composition of bacterial components
bcomps <- components(grafs[["Bacteria"]])
bmems <- membership(bcomps)
grafs[["Bacteria"]] %>% 
  activate(nodes) %>% 
  mutate(component = as.character(bmems[name])) %>% 
  create_layout("manual", node.positions = as.data.frame(select(., x, y))) %>% 
  plot.mapper(aes_(size = ~size, color = ~component))
save_plot(paste0(figs.dir, "nahant-bac-components.pdf"), last_plot(),
          base_height = 6)
v2p[["Bacteria"]][, component := bmems[vertex.name]]
setkey(nahant, kingdom)
phyla.comps <- merge(nahant["Bacteria"], v2p[["Bacteria"]], by = "day",
                     allow.cartesian = TRUE) %>% 
  .[, .(mean.freq = mean(freq)), by = .(component, phylum)]
phyla.comps[, rank := frank(-mean.freq), by = component]

#' Let's worry only about the 2 largest components:
setkey(phyla.comps, component)
phyla.comps[.(c(1, 2))] %>% 
  dcast(phylum ~ component, value.var = "mean.freq") %>% 
  ggplot(aes(x = `1`, y = `2`)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  geom_point(aes(color = phylum), 
             data = function(d) filter(d, `1` > 0.00001 & `2` > 0.000025)) +
  scale_color_brewer(palette = "Dark2")
save_plot(paste0(figs.dir, "nahant-bac-components-phyla.pdf"), last_plot(),
          base_height = 6)

#' # Diatom vs dinoflagellate dominance
setkey(nahant, order)
dd <- nahant[c("Diatomea", "Dinoflagellata"), .(freq = sum(freq)),
             by = .(order, day)]
dd <- dcast(dd, day ~ order, value.var = "freq")
dd[, ratio := Diatomea / Dinoflagellata]
hist(log10(dd$ratio))
v2p <- lapply(v2p, merge, y = dd, by = "day")
vdd <- lapply(v2p, function(d) {
  d[, .(mean.dia = mean(Diatomea),
        mean.dino = mean(Dinoflagellata),
        mean.ratio = mean(ratio)), by = vertex]
})
grafs <- mapply(left_join, x = grafs, y = vdd,
                MoreArgs = list(by = c("vertex" = "vertex")), SIMPLIFY = FALSE)
grafs <- lapply(grafs, function(g) {
  g %>%
    activate(nodes) %>% 
    mutate(euk.type = sapply(mean.ratio, function(r) {
      lr <- log10(r)
      if (lr > 0) {
        "dia"
      } else {
        "dino"
      }
    }))
})
dd.plots <- lapply(grafs, function(g) {
  g <- mutate(g, log.mean.ratio = log10(mean.ratio))
  lo <- create_layout(g, "manual",
                      node.positions = as.data.frame(select(g, x, y)))
  # plot.mapper(lo, aes_(size = ~size, color = ~log.mean.ratio)) +
  #   scale_color_gradient2(low = "blue", high = "red") +
  #   labs(color = "log10(Dia/Dino)")
  plot.mapper(lo, aes_(size = ~size, color = ~euk.type)) +
    labs(color = "dominant")
})
plot_grid(plotlist = dd.plots, labels = names(dd.plots))
# save_plot(paste0(figs.dir, "nahant-log-dd-ratio.pdf"), last_plot(),
#           ncol = 2, base_height = 6)
save_plot(paste0(figs.dir, "nahant-dominant-euk.pdf"), last_plot(),
          ncol = 2, base_height = 6)


# frames ------------------------------------------------------------------

grafs <- mapply(function(graf, lo) {
  graf %>%
    activate(nodes) %>%
    mutate(x = lo$x, y = lo$y)
}, graf = grafs, lo = layouts, SIMPLIFY = FALSE)
sample.subgraphs <- mapply(function(graf, v2p) {
  spl <- split(v2p, by = "day")
  lapply(spl, function(s, graf) {
    graf %>%
      slice(s$vertex) %>%
      mutate(day = s$day)
  }, graf = graf)
}, graf = grafs, v2p = v2p, SIMPLIFY = FALSE)
sample.plots <- mapply(function(sgs, layout) {
  lapply(sgs, function(sg, layout) {
    sg <- as.data.frame(sg)
    day <- unique(sg$day)
    p <- plot.mapper(layout, aes_(size = ~size), NULL, color = "grey")
    p + geom_point(aes(x = x, y = y, size = size),
                 data = sg, color = "blue") +
      labs(title = day)
  }, layout = layout)
}, sgs = sample.subgraphs, layout = layouts)
for (dt in samples) {
  setorder(dt, day)
  dt[, id := seq_len(.N)]
}
for (kingdom in names(samples)) {
  for (x in seq_len(nrow(samples[[1]]))) {
    day <- samples[[kingdom]][x, day]
    fn <- paste0(kingdom, x, ".png")
    save_plot(paste0("nahant-trajectories/frames/", fn),
              sample.plots[[kingdom]][[x]])
  }
}

#' # Joint bac-euk analysis
# jdays <- sort(rownames(distances[["Bacteria"]]))
# jdist <- distances %>%
#   lapply(function(m, jdays) m[jdays, jdays], jdays = jdays) %>%
#   lapply(function(x) x ^ 2) %>%
#   do.call("+", .) %>%
#   sqrt
# jmds2 <- cmdscale(jdist, eig = TRUE)
# jmds2$GOF
# jmds2$points %>%
#   as.data.table(keep.rownames = "day") %>%
#   .[, day := as.numeric(day)] %>%
#   ggplot(aes(x = V1, y = V2)) +
#   geom_point(aes(color = day)) +
#   scale_color_distiller(palette = "Spectral")
# jrkmds2 <- jmds2$points %>% as.data.frame %>% lapply(rank)

# joint mapper ------------------------------------------------------------

# jmpr <- mapper2D(jdist, jrkmds2, num_intervals = ni, percent_overlap = po)
# jv2p <- vertex.2.points(jmpr$points_in_vertex)
# ni <- c(10, 10)
# po <- 70
# jmpr <- mapper2D(jdist, jrkmds2, num_intervals = ni, percent_overlap = po)
# jv2p <- vertex.2.points(jmpr$points_in_vertex)
# jv2p[, day := as.numeric(jdays[point])]
# jverts <- jv2p[, .(size = .N, mean.day = mean(day)),
#                by = .(vertex, vertex.name)]
# jgraf <- mapper.2.igraph(jmpr)
# jgraf <- jgraf %>%
#   as_tbl_graph %>%
#   activate(nodes) %>%
#   left_join(jverts, by = c("name" = "vertex.name"))
# set.seed(0)
# lo <- create_layout(jgraf, "fr")
# jgraf <- mutate(jgraf, x = lo$x, y = lo$y)
# plot.mapper(lo, aes_(size = ~size, color = ~mean.day)) +
#   scale_color_distiller(palette = "Spectral")

#' ## Variance of phyla across time
# phyla.days <- nahant[, .(freq = sum(freq)), by = .(kingdom, phylum, day)]
# phyla <- phyla.days[, .(mean = mean(freq), variance = var(freq)),
#                     by = .(kingdom, phylum)]

#' Variance scales with mean abundance:
# plot(phyla$mean, phyla$variance)

#' But the most variable phyla are not the rarest so that's good (not noisy).
# phyla[, scaled.var := variance / mean]
# phyla[, scaled.rk := frank(-scaled.var), by = kingdom]
# ggplot(phyla, aes(x = scaled.rk, y = mean)) +
#   geom_point(data = function(d) filter(d, scaled.rk > 5)) +
#   geom_point(aes(color = phylum),
#              data = function(d) filter(d, scaled.rk <= 5)) +
#   scale_y_log10() +
#   facet_wrap(~ kingdom, nrow = 2)
#
# phyla.verts <- merge(phyla.days, jv2p, by = "day", allow.cartesian = TRUE) %>%
#   .[, .(mean = mean(freq)), by = .(kingdom, phylum, vertex)]
# phyla.verts <- phyla.verts %>%
#   filter(!is.na(phylum)) %>%
#   dcast(vertex ~ phylum, value.var = "mean")
# jgraf <- jgraf %>%
#   activate(nodes) %>%
#   left_join(phyla.verts, by = c("vertex" = "vertex"))
# phyla.plots <- phyla[scaled.rk <= 5] %>%
#   split(by = "kingdom") %>%
#   lapply(function(k, jgraf) {
#     ps <- k[!is.na(phylum), phylum]
#     lapply(ps, function(ph, jgraf) {
#       xy <- as.data.frame(jgraf)[, c("x", "y")]
#       lo <- create_layout(jgraf, "manual", node.positions = xy)
#       plot.mapper(lo, aes_string(size = "size", color = ph)) +
#         scale_color_distiller(palette = "Spectral")
#     }, jgraf = jgraf)
#   }, jgraf = jgraf)

#' ## Mean abundance of most variable bacterial phyla across the joint space
# phyla.plots[["Bacteria"]]

#' ## Mean abundance of most variable eukaryotic phyla across the joint space
# phyla.plots[["Eukaryota"]]
