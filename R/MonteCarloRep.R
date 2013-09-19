MonteCarloRep <-
function (x.matrix, ind, nit = 100, ComparisonFunc = NULL)
  # Calculates x.matrix repeatability using parametric sampling
  #
  # Args:
  #     x.matrix: covariance or correlation matrix.
  #               if x.matrix is a correlation matrix will use MantelCor,
  #               else, will use RandomSkewers
  #     ind: number of indivuals on each sample
  #     nit: number of samples
  #     Comparisonfunc: Arbitrary function to compare 2 matrices. Must return single numeric value.
  # Return:
  #     mean correlation of sample covariance matrices with original input x.matrix
{
  if(!require(mvtnorm)) install.packages("mvtnorm")
  library(mvtnorm)
  if (sum(diag(x.matrix)) == dim (x.matrix) [1]){
    Func <- function(x, y, z) MantelCor(x, y, z)[1]
    Type <- cor
  }
  else{
    Func <- function(x, y, z) RandomSkewers(x, y, z)[1]
    Type <- var
  }
  if(!is.null(ComparisonFunc)) Func <- function(x, y, z) ComparisonFunc(x, y)
  R <- c()
  for (N in 1:nit){
    rand.samp <- rmvnorm (ind, rep(0, times = dim (x.matrix)[1]),
                           sigma = x.matrix, method = "chol")
    rand.matrix <- Type (rand.samp)
    R[N] <- Func (x.matrix, rand.matrix, 1000)
  }
  return (mean(R))
}
