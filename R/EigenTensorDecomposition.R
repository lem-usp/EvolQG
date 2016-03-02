#' Eigentensor Decomposition 
#' 
#' This function performs eigentensor decomposition on a set of covariance matrices.
#'
#' @param matrix.array k x k x m array of m covariance matrices with k traits;
#' @param return.projection Should we project covariance matrices into estimated eigentensors? Defaults to TRUE.
#' @return List with the following components:
#' @return values vector of ordered eigenvalues associated with eigentensors;
#' @return matrices array of eigentensor in matrix form;
#' @return projection matrix of scaled projected covariance matrices over eigentensors.
#' 
#' @details The number of estimated eigentensors is the minimum between k(k+1)/2 - 1 and m - 1, 
#' in a similar manner to the usual principal component analysis.
#' 
#' @references Basser P. J., Pajevic S. 2007. Spectral decomposition of a 4th-order 
#' covariance tensor: Applications to diffusion tensor MRI. Signal Processing. 87:220–236.
#' @references Hine E., Chenoweth S. F., Rundle H. D., Blows M. W. 2009. Characterizing 
#' the evolution of genetic variance using genetic covariance tensors. Philosophical 
#' transactions of the Royal Society of London. Series B, Biological sciences. 364:1567–78.
#' 
#' @author Guilherme Garcia
#' @seealso \code{\link{ProjectMatrix}}, \code{\link{RevertMatrix}}
#' 
#' @examples 
#' data(dentus)
#'
#' dentus.vcv <- daply (dentus, .(species), var) [, -5, -5]
#' # There will be some warnings
#'
#' dentus.vcv <- aperm(dentus.vcv, c(2, 3, 1))
#'
#' dentus.etd <- EigenTensorDecomposition(dentus.vcv, TRUE)
#'
#' # Plot some results
#' par(mfrow = c(1, 2))
#' plot(dentus.etd $ values, pch = 20, type = 'b', ylab = 'Eigenvalue')
#' plot(dentus.etd $ projection [, 1:2], pch = 20)
#' text(dentus.etd $ projection [, 1:2],
#'      labels = rownames (dentus.etd $ projection), pos = 2)
#' 
#' @importFrom matrixcalc frobenius.prod
#' @importFrom expm logm sqrtm
#' 
#' @export
#'
EigenTensorDecomposition <-
  function (matrix.array, return.projection = TRUE)
  {
    mean.matrix <- MeanMatrix(matrix.array)
    mean.is <- sqrtm(solve(mean.matrix))
    lce.array <- aaply(matrix.array, 3, 
                       function (A) logm(mean.is %*% A %*% mean.is))  
    lce.array <- aperm(lce.array, c(2, 3, 1))
    n.traits <- dim(lce.array) [1]
    n.matrix <- dim(lce.array) [3]
    
    ### Build Sigma
    variances <- aaply (matrix.array, 3, diag)
    covariances <- aaply (matrix.array, 3, function (x) x [lower.tri (x)])
    block.var <- var (variances)
    block.off <- cov (variances, covariances)
    block.cov <- var (covariances)
    upper.block <- cbind (block.var, sqrt (2) * block.off)
    lower.block <- cbind (sqrt (2) * t (block.off), 2 * block.cov)
    Sigma <- rbind (upper.block, lower.block)
    
    n.tensor <- min (n.matrix - 1, dim (Sigma) [1])
    eigen.dec <- eigen (Sigma)
    eigen.matrices <-
      aaply (eigen.dec $ vectors [, 1:n.tensor], 2, ### ematrices with non-zero evals
             function (x)
             {
               eigen.mat <- diag (x [1:n.traits])
               eigen.mat [upper.tri (eigen.mat)] <-
                 x [(n.traits+1):length (x)] / sqrt (2)
               eigen.mat <- eigen.mat + t(eigen.mat)
               diag (eigen.mat) <- diag (eigen.mat) / 2
               eigen.mat
             })
    eigen.matrices <- aperm (eigen.matrices, c(2, 3, 1), resize = TRUE)
    dimnames (eigen.matrices) <-
      list (rownames (matrix.array),
            colnames (matrix.array),
            paste ('Eigentensor', 1:n.tensor, sep = ' '))
    
    out <- list ('mean' = mean.matrix,
                 'values' = eigen.dec $ values [1:n.tensor],
                 'matrices' = eigen.matrices)
   
    if (return.projection)
    {
      project <- aaply (lce.array, 3, 
                        function (A, B) aaply (B, 3, frobenius.prod, y = A),
        B = eigen.matrices)
      out $ projection
      ### already centered because of matrix product with the mean
      out $ projection <- scale (out $ projection, center = FALSE, scale = TRUE)
      ## out $ projection <- data.frame(out $ projection)
      ## if(!is.null(dimnames(matrix.array) [3]))
      ##  out $ projection <- data.frame('otu' = dimnames(matrix.array) [3], out $ projection)
    }
    return (out)
  }
