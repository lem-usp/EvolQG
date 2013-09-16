RandomSkewers <-
function (cov.matrix.1, cov.matrix.2, nsk = 10000)
  # Calculates covariance matrix correlation via random skewers
  # Args:
  #     cov.matrix.(1,2): Two covariance matrices to be compared
  #     nsk: Number of generated random skewers
  # Return:
  #     List with mean value of correlation, p value and standard deviation
{
  traits <- dim (cov.matrix.1) [1]
  base.vector <- Normalize(rnorm(traits))
  random.vectors <- array (rnorm (nsk * traits, mean = 0, sd = 1), c(traits, nsk))
  random.vectors <- apply (random.vectors, 2, Normalize)
  dist <- base.vector %*% random.vectors
  dz1 <- apply (cov.matrix.1 %*% random.vectors, 2, Normalize)
  dz2 <- apply (cov.matrix.2 %*% random.vectors, 2, Normalize)
  real <- apply (dz1 * dz2, 2, sum)
  ac <- mean (real)
  stdev <- sd (real)
  prob <- sum (ac < dist) / nsk
  output <- c(ac, prob, stdev)
  names(output) <- c("AC","P","SD")
  return(output)
}
