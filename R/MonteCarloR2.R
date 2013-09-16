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
  #   a vector with iterations + 1 entries; the first entry covesponds
  #   to the actual r2 value calculated from the original matrix. The
  #   remaining entries are values calculated from the resampling procedure.
{
  n.traits <- dim (cov.matrix) [1]
  populations <- list ()
  for (i in 1:iterations)
    populations [[i]] <- rmvNorm2 (sample.size, sigma = cov.matrix, method = 'chol')
  it.matrices <- lapply (populations, cor)
  it.r2 <- sapply (it.matrices, CalcR2)
  r2 <- c (CalcR2 (cov2cor(cov.matrix)), it.r2)
  return (r2)
}
