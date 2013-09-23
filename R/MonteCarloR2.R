MonteCarloR2 <-
function (cov.matrix, sample.size, iterations = 1000)
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
  n.traits <- dim (cov.matrix) [1]
  populations <- list ()
  for (i in 1:iterations)
    populations [[i]] <- rmvnorm (sample.size, sigma = cov.matrix, method = 'chol')
  it.matrices <- lapply (populations, cor)
  it.r2 <- sapply (it.matrices, CalcR2)
  return (it.r2)
}
