MultiMahalanobis <- function (mean.list, cov.matrix, num.cores = 1)
  # Calculates the Mahalanobis distance between a list of species mean,
  # using a global covariance matrix
  #
  # Args:
  #  mean.list: a list of species mean
  #  cov.matrix: global covariance matrix
  #
  # Return:
  #  returns a matrix of species-species comparisons.
{
  library(plyr)
  library(reshape2)
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else{
    parallel = FALSE
  }
  num.mean <- length(mean.list)
  if(is.null(names(mean.list))) {names(mean.list) <- 1:num.mean}
  mean.names <- names (mean.list)
  CompareToN <- function(n) ldply(mean.list[(n+1):num.mean],
                                  function(x) {mahalanobis(x, mean.list[[n]], cov.matrix)},
                                  .parallel = parallel)
  comparisons <- adply(1:(num.mean-1), 1,  CompareToN, .parallel = parallel)
  dists <- acast(comparisons[-4], X1~.id)[,mean.names[-1]]
  distances <- array (0, c(num.mean, num.mean))
  distances[upper.tri(distances)] <- dists[upper.tri(dists, diag=T)]
  rownames (distances) <- mean.names
  colnames (distances) <- mean.names
  return (as.dist(distances))
}
