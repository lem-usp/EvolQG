MonteCarloR2 <- function (cov.matrix, sample.size, iterations = 1000, num.cores = 1)
  # Computes a distribution of magnitudes of integration (r2)
  # for a given covariance matrix.
  #
  # Args:
  #   cov.matrix: a square symmetric covariance matrix
  #   sample.size: number of individuals to sample
  #   iterations: number of populations sampled
  # Return:
  # a vector whose entries are R2 values calculated from the
  # resampling procedure.
{
  library(mvtnorm)
  library(plyr)
  library(reshape2)
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else
    parallel = FALSE
  n.traits <- dim (cov.matrix) [1]
  populations  <- alply(1:iterations, 1,
                        function(x) rmvnorm (sample.size, sigma = cov.matrix, method = 'chol'),
                        .parallel=parallel)
  it.r2 <- laply (populations, function (x) CalcR2(cor(x)), .parallel = parallel)
  return (it.r2)
}
