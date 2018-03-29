library(data.table)
library(TDA)
library(philentropy)

util.dir <- "../r/"

source(paste0(util.dir, "load_cholera_data.R"))

gordon[, freq := count / sum(count), by = sample]
distribs <- dcast(gordon, sample ~ otu, value.var = "freq", fill = 0)
sample.names <- distribs$sample
distribs <- as.matrix(distribs[, -1])
jsd <- JSD(distribs)
rownames(jsd) <- sample.names
colnames(jsd) <- sample.names
js.dist <- sqrt(jsd)

# density estimator on a grid method
# mds
n <- length(rn)
k <- 2
x <- cmdscale(js.dist, k = k, eig = TRUE) 
print("MDS goodness of fit:")
print(x$GOF)
mds <- x$points
ngrid <- 100
grid <- expand.grid(
  as.data.frame(
  apply(mds, 2, function(v) seq(min(v), max(v), length.out = ngrid)
)))

#maxscale <- max(js.dist)
#print(paste("Maxscale = max distance =" , maxscale))
#
#rips <- ripsDiag(js.dist, maxdimension = 0, maxscale = maxscale,
#  dist = "arbitrary", library = "Dionysus")

# pdf("figs/cholera_persistence_diag.pdf")
# plot(rips$diagram, main = "Diagram")
# plot(rips$diagram, barcode = TRUE, main = "Barcode")
# cdf <- ecdf(rips$diagram[, 3])
# plot(cdf, main = "b0")
# bzero <- function(x) (1 - cdf(x)) * nrow(rips$diagram)
# plot(bzero, main = "# components", log = "y")
# dev.off()
