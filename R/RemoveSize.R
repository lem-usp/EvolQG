RemoveSize <- function (cov.matrix) {
  cov.matrix.svd  <-  svd(cov.matrix)
  size.eigen.vector <- cov.matrix.svd$u[, 1]
  size.eigen.value <- cov.matrix.svd$d[1]
  size.factor <- size.eigen.vector * sqrt(size.eigen.value)
  cov.matrix.size.removed <- cov.matrix - size.factor %*% t(size.factor)
  return (cov.matrix.size.removed)
}
