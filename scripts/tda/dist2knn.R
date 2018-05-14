dist2knn <- function(dist.mat, k) {
  source("k.first.R")
  knn <- apply(dist.mat, 1, k.first, k = k)
  names(knn) <- rownames(dist.mat)
  knn
}