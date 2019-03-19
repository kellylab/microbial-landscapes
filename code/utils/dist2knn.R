dist2knn <- function(dist.mat, k) {
  k.first <- function(v, k, decreasing = FALSE) {
    v <- sort(v[v > 0], decreasing = decreasing)
    mean(v[1:k])
  }
  knn <- apply(dist.mat, 1, k.first, k = k)
  names(knn) <- rownames(dist.mat)
  knn
}