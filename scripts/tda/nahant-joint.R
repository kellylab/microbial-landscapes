library(tidyverse)
library(data.table)
library(TDAmapper)
library(ggraph)
library(igraph)
library(tidygraph)
library(cowplot)

# setup -------------------------------------------------------------------


scripts.dir <- "../r/"
source(paste0(scripts.dir, "load-nahant-data.R"))
source("k.first.R")
source("vertex.2.points.R")
source("mapper.2.igraph.R")
source("plot.mapper.R")
source("dist2knn.R")

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

#' # Joint bac-euk analysis
jdays <- sort(rownames(distances[["Bacteria"]]))
jdist <- distances %>%
  lapply(function(m, jdays) m[jdays, jdays], jdays = jdays) %>%
  lapply(function(x) x ^ 2) %>%
  do.call("+", .) %>%
  sqrt
k <- 10
knn <- dist2knn(jdist, k)
hist(knn)
jmds2 <- cmdscale(jdist, eig = TRUE)
jmds2$GOF
jrkmds2 <- apply(jmds2$points, 2, rank)
jrkmds2 %>%
  as.data.table(keep.rownames = "day") %>%
  .[, day := as.numeric(day)] %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(color = day)) +
  scale_color_distiller(palette = "Spectral")

# joint mapper ------------------------------------------------------------

ni <- c(10, 10)
po <- 70
jmpr <- mapper2D(jdist, as.data.frame(jrkmds2),  num_intervals = ni,
                 percent_overlap = po)
jv2p <- vertex.2.points(jmpr$points_in_vertex)
jv2p[, day := as.numeric(jdays[point])]
jv2p[, knn := knn[as.character(day)]]
jverts <- jv2p[, .(size = .N, mean.day = mean(day), mean.knn = mean(knn)),
               by = .(vertex, vertex.name)]
jgraf <- mapper.2.igraph(jmpr)
jgraf <- jgraf %>%
  as_tbl_graph %>%
  activate(nodes) %>%
  left_join(jverts, by = c("name" = "vertex.name"))
set.seed(0)
lo <- create_layout(jgraf, "fr")
jgraf <- mutate(jgraf, x = lo$x, y = lo$y)
plot.mapper(lo, aes_(size = ~size, color = ~mean.day)) +
  scale_color_distiller(palette = "Spectral")
plot.mapper(lo, aes_(size = ~size, color = ~mean.knn)) +
  scale_color_distiller(palette = "Blues")

#' ## Variance of phyla across time
phyla.days <- nahant[, .(freq = sum(freq)), by = .(kingdom, phylum, day)]
phyla <- phyla.days[, .(mean = mean(freq), variance = var(freq)),
                    by = .(kingdom, phylum)]

#' Variance scales with mean abundance:
plot(phyla$mean, phyla$variance)

#' But the most variable phyla are not the rarest so that's good (not noisy).
phyla[, scaled.var := variance / mean]
phyla[, scaled.rk := frank(-scaled.var), by = kingdom]
ggplot(phyla, aes(x = scaled.rk, y = mean)) +
  geom_point(data = function(d) filter(d, scaled.rk > 5)) +
  geom_point(aes(color = phylum),
             data = function(d) filter(d, scaled.rk <= 5)) +
  scale_y_log10() +
  facet_wrap(~ kingdom, nrow = 2)

#' ## Variance of phyla across joint space
phyla.verts <- merge(phyla.days, jv2p, by = "day", allow.cartesian = TRUE) %>%
  .[, .(mean = mean(freq)), by = .(kingdom, phylum, vertex)]
phyla.verts <- phyla.verts %>%
  filter(!is.na(phylum)) %>%
  dcast(vertex ~ phylum, value.var = "mean")
jgraf <- jgraf %>%
  activate(nodes) %>%
  left_join(phyla.verts, by = c("vertex" = "vertex"))
phyla.plots <- phyla[scaled.rk <= 5] %>%
  split(by = "kingdom") %>%
  lapply(function(k, jgraf) {
    ps <- k[!is.na(phylum), phylum]
    lapply(ps, function(ph, jgraf) {
      xy <- as.data.frame(jgraf)[, c("x", "y")]
      lo <- create_layout(jgraf, "manual", node.positions = xy)
      plot.mapper(lo, aes_string(size = "size", color = ph)) +
        scale_color_distiller(palette = "Spectral")
    }, jgraf = jgraf)
  }, jgraf = jgraf)

#' ### Mean abundance of most variable bacterial phyla across the joint space
phyla.plots[["Bacteria"]]

#' ### Mean abundance of most variable eukaryotic phyla across the joint space
phyla.plots[["Eukaryota"]]

#' ## Covariance of phyla across the joint space
phyla.cor <- phyla.verts[, -1] %>%
  as.matrix %>%
  cor

#' Inverting the Pearson correlation matrix gives coupling matrix constrained by
#' the correlations while maximizing entropy
inv <- solve(phyla.cor)
hist(inv)
setkey(phyla, kingdom)
bphy <- phyla["Bacteria", phylum]
bphy <- bphy[!is.na(bphy)]
ephy <- phyla["Eukaryota", phylum]
ephy <- ephy[!is.na(ephy)]
inv[bphy, ephy] %>%
  melt(varnames = c("bacteria", "eukaryota"), value.name = "coupling") %>%
  ggplot(aes(x = bacteria, y = eukaryota)) +
  geom_tile(aes(fill = coupling)) +
  coord_equal() +
  scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
setkey(phyla, phylum)
coupling.graf <- graph_from_adjacency_matrix(inv, mode = "upper",
                                             weighted = TRUE, diag = FALSE) %>%
  as_tbl_graph %>%
  activate(nodes) %>%
  mutate(kingdom = phyla[name, kingdom]) %>%
  activate(edges) %>%
  mutate(positive = weight > 0, weight = abs(weight))
ggraph(coupling.graf, "fr", niter = 1000) +
  geom_edge_link0(aes(color = positive, alpha = weight)) +
  geom_node_point(aes(color = kingdom)) +
  theme_graph()
