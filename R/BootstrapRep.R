BootstrapRep <-
function (ind.data, nb = 100)
  # Calculates the repeatability of the covariance matrix of the suplied data
  # via bootstrap ressampling
  #
  # Args:
  #     ind.data: original individual data
  #     nb = number of resamples
  # Return:
  #     returns the mean repeatability
{
  n.ind <-  dim (ind.data) [1]
  original.cov.matrix <- var (ind.data)
  v.rep <- c()
  for (N in 1:nb){
    sampled.data <- sample (1:n.ind, n.ind, TRUE)
    sampled.data.cov.matrix <- var (ind.data[sampled.data,])
    v.rep [N] <- RandomSkewers (original.cov.matrix, sampled.data.cov.matrix, 1000) [1]
  }
  out <- mean (v.rep)
  return (out)
}
