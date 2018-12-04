height.distribs <- function(distances, points.in.level, nbin) {
  get.lvldist <- function(pil, distances) {
    as.dist(distances[pil, pil])
  }
  get.lengths <- function(heights, distances) {
    c(heights, max(distances))
  }
  bin <- function(lengths, nbin) {
    breaks <- seq(from = min(lengths), to = max(lengths),
                  by = (max(lengths) - min(lengths)) / nbin)
    h <- hist(lengths, breaks = breaks, plot = FALSE)
    data.table(bin = seq_along(h$counts), n = h$counts)
  }
  lvldists <- lapply(points.in.level, get.lvldist, distances = distances)
  clusts <- lapply(lvldists, hclust)
  heights <- lapply(clusts, function(c) c$height)
  lengths <- mapply(get.lengths, heights = heights, distances = lvldists)
  bins <- lapply(heights, bin, nbin = nbin)
  rbindlist(bins, idcol = "level")
}