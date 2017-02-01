#' Mean Covariance Matrix
#'
#' Estimate geometric mean for a set of covariance matrices
#'
#' @param matrix.array k x k x m array of covariance matrices, with k traits and m matrices
#' @param tol minimum riemannian distance between sequential iterated means for accepting an estimated matrix
#' @return geometric mean covariance matrix
#' 
#' @importFrom matrixcalc frobenius.norm
#' @importFrom Matrix Schur
#' 
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
MeanMatrix <- function (matrix.array, tol = 1e-10)
{
  A <- matrix.array
  m <- dim(A) [3]

  niold <- Inf
  
  A.chol <- aaply(A, 3, chol)
  X <- aaply (A, c(1, 2), mean)
  
  repeat
  {
    X.chol <- chol(X)
    X.chinv <- solve(X.chol)
    UV <- alply(A.chol, 1, function(R) 
    {
      Z <- R %*% X.chinv
      Schur(t(Z) %*% Z, vectors = TRUE)
    })
    V <- laply(UV, function(L) L $ EValues)
    U <- laply(UV, function(L) L $ Q)
    
    ## dynamic choice of theta
    ch = aaply(V, 1, max) / aaply(V, 1, min)
    dh = log(ch) / (ch - 1)
    beta = sum(ch)
    gamma = sum(ch * dh)
    theta = 2 / (gamma + beta)
    
    Ti <- aaply(1:m, 1, function(i)
      {
        T = U[i, , ] %*% diag(log(V[i, ])) %*% t(U [i, , ])
        T + t(T) / 2
      })
    S = aaply(Ti, c(2, 3), sum)
    UV.S <- Schur(S)
    Z <- diag(exp(UV.S $ EValues * theta / 2)) %*% t(UV.S $ Q) %*% X.chol
    
    X.next <-t(Z) %*% Z
    
    ni <- max(abs(UV.S $ EValues))
    
    if(ni < frobenius.norm(X.next) * tol || ni > niold)
      break
    else
    {
      niold <- ni 
      X <- X.next
    }

  }
  X.next
}
