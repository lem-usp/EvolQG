#' Mean Covariance Matrix
#'
#' Estimate geometric mean for a set of covariance matrices
#'
#' @param matrix.array k x k x m array of covariance matrices, with k traits and m matrices
#' @param tol minimum riemannian distance between sequential iterated means for accepting an estimated matrix
#' @param verbose print values for each iteration?
#' @return geometric mean covariance matrix
#' 
#' @importFrom expm sqrtm logm expm
#' @rdname MeanMatrix
#' @references Bini, D. A., Iannazzo, B. 2013. Computing the Karcher Mean of Symmetric
#' Positive Definite Matrices. Linear Algebra and Its Applications, 16th ILAS Conference
#' Proceedings, Pisa 2010, 438 (4): 1700-1710. doi:10.1016/j.laa.2011.08.052. 
#' @author Guilherme Garcia, Diogo Melo
#' 
#' @seealso \code{\link{EigenTensorDecomposition}}, \code{\link{RiemannDist}}
#'
#' @export
#' 
MeanMatrix <- function (matrix.array, tol = 1e-10, verbose = FALSE)
{
  A <- matrix.array
  m <- dim(A) [3]
  
  v <- 1/m
  
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
    
    if(verbose) print(dF)
    if(dF < tol) break
    
    X <- X.next
  }
  X.next
}
