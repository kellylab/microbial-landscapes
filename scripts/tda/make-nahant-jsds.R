library(data.table)
library(philentropy)

scripts.dir <- "../r/"

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