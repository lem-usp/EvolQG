RemoveSize <-
function (cov.matrix)
  # Removes first principal component effect in cov.matrix.
  #
  # Args:
  #   cov.matrix: A simetric covariance matrix
  # Return:
  #   cov.matrix with size removed
{
  cov.matrix.svd  <-  svd(cov.matrix)
  size.eigen.vector <- cov.matrix.svd$u[, 1]
  size.eigen.value <- cov.matrix.svd$d[1]
  size.factor <- size.eigen.vector * sqrt(size.eigen.value)
  cov.matrix.size.removed <- cov.matrix - size.factor %*% t(size.factor)
  return (cov.matrix.size.removed)
}
