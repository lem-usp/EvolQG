KrzCor <-
function (cov.matrix.1, cov.matrix.2, ret.dim = NULL)
  # Calculates the Krzanowski correlation between matrices
  #
  # Args:
  #     cov.matrix.(1,2): covariance being compared
  #     ret.dim: number of retained dimensions in the comparison,
  #              default for nxn matrix is n/2-1
  # Return:
  #     Kzranowski correlation
{
  if (is.null(ret.dim))
    ret.dim = dim(cov.matrix.1)[1]/2 - 1
  EigenVectors <- function (x) return (eigen(x)$vectors[,1:ret.dim])
  A <- EigenVectors (cov.matrix.1)
  B <- EigenVectors (cov.matrix.2)
  S <- t(A) %*% B %*% t(B) %*% A
  SL <- sum (eigen(S)$values) / ret.dim
  return (SL)
}
