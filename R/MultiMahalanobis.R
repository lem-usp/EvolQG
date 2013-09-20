MultiMahalanobis <-
function (mean.list, cov.matrix)
  # Calculates the Mahalanobis distance between a list of species means,
  # using a global covariance matrix
  #
  # Args:
  #  mean.list: a list of species means
  #  cov.matrix: global covariance matrix
  #
  # Return:
  #  returns a matrix of species-species comparisons.
{
  n.means <- length (mean.list)
  means.names <- names (mean.list)
  distances <- array (0, c(n.means, n.means))
  for (i in 1:(n.means - 1)) {
    for (j in (i+1):n.means) {
      comparing.now <- mahalanobis (mean.list [[i]],
                                    mean.list [[j]],
                                    cov.matrix)
      distances [i, j] <- distances [j, i] <- comparing.now
    }
  }
  rownames (distances) <- means.names
  colnames (distances) <- means.names
  return (as.dist(distances))
}
