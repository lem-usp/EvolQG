#' Covariance Tensor between Covariance Matrices
#'
#' Computes the two-dimensional projection of the fourth-order 
#' covariance tensor between covariance matrices.
#' 
#' @param matrix.array m x m x k array of k covariance matrices with m traits
#' @return projected covariance tensor: a square symmetric matrix with m(m + 1)/2 rows
#' @rdname BuildSigma
#' @references Basser P. J., Pajevic S. 2007. Spectral decomposition of a 4th-order 
#' covariance tensor: Applications to diffusion tensor MRI. Signal Processing. 87:220–236.
#' @references Hine E., Chenoweth S. F., Rundle H. D., Blows M. W. 2009. Characterizing 
#' the evolution of genetic variance using genetic covariance tensors. Philosophical 
#' transactions of the Royal Society of London. Series B, Biological sciences. 364:1567–78.
#'
#' @author Guilherme Garcia
#' @seealso \code{\link{EigenTensorDecomposition}}
#' 
BuildSigma <- function (matrix.array)
{
  variances <- aaply (matrix.array, 3, diag)
  covariances <- aaply (matrix.array, 3, function (x) x [lower.tri (x)])
  block.var <- var (variances)
  block.off <- cov (variances, covariances)
  block.cov <- var (covariances)
  upper.block <- cbind (block.var, sqrt (2) * block.off)
  lower.block <- cbind (sqrt (2) * t (block.off), 2 * block.cov)
  Sigma <- rbind (upper.block, lower.block)
  return (Sigma)
}
