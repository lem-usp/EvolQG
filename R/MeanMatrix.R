#' Mean Covariance Matrix
#'
#' Estimate geometric mean for a set of covariance matrices
#'
#' @param matrix.array m x m x k array of covariance matrices, with m traits and k matrices
#' @param tol minimum value between iterations for evaluating convergence
#' @return geometric mean covariance matrix
#' 
#' @importFrom expm sqrtm logm expm
#' @rdname MeanMatrix
#' @references Bini, D. A., Iannazzo, B. 2013. Computing the Karcher Mean of Symmetric
#' Positive Definite Matrices. Linear Algebra and Its Applications, 16th ILAS Conference
#' Proceedings, Pisa 2010, 438 (4): 1700–1710. doi:10.1016/j.laa.2011.08.052. 
#' @author Guilherme Garcia
#' @seealso \code{\link{EigenTensorDecomposition}}
#'
MeanMatrix <- function (matrix.array, tol = 1e-10)
{
  A <- matrix.array
  k <- dim(A) [3]
  
  v <- 1/k
  
  ## from here, each matrix is the first dimension slice of the array
  A.inv <- aaply(A, 3, solve)
  
  X <- aaply (A, c(1, 2), mean)
  
  repeat
  {
    
    X.sq <- sqrtm (X)
    
    A.prod <- aaply (A.inv, 1, function (Ai) logm(X.sq %*% Ai %*% X.sq))
    
    A.sum <- aaply (A.prod, c(2, 3), sum)
    
    X.next <- X.sq %*% expm(- v * A.sum) %*% X.sq
    
    dF <- RiemannDist(X, X.next)
    if(dF < tol)
      break
    
    print (dF)
    X <- X.next
  }
  X.next
}
