#' Eigentensor Decomposition 
#' 
#' This function performs eigentensor decomposition on a set of covariance matrices
#'
#' @param matrix.array k x k x m array of m covariance matrices with k traits;
#' @param return.projection Should we project covariance matrices into estimated eigentensors? Defaults to TRUE.
#' @return List with the following components:
#' @return values vector of ordered eigenvalues associated with eigentensors;
#' @return matrices array of eigentensor in matrix form;
#' @return projection matrix of projected covariance matrices over eigentensors.
#' 
#' @details The number of estimated eigentensors is the minimum between k(k+1)/2 and m, 
#' in a similar manner to the usual principal component analysis.
#' 
#' @references Basser P. J., Pajevic S. 2007. Spectral decomposition of a 4th-order 
#' covariance tensor: Applications to diffusion tensor MRI. Signal Processing. 87:220–236.
#' @references Hine E., Chenoweth S. F., Rundle H. D., Blows M. W. 2009. Characterizing 
#' the evolution of genetic variance using genetic covariance tensors. Philosophical 
#' transactions of the Royal Society of London. Series B, Biological sciences. 364:1567–78.
#' 
#' @importFrom matrixcalc frobenius.prod 
#' @export

EigenTensorDecomposition <-
  function (matrix.array, return.projection = FALSE)
  {
    n.traits <- dim (matrix.array) [1]
    n.matrix <- dim (matrix.array) [3]
    Sigma <- BuildSigma (matrix.array)
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
            paste ('PM', 1:n.tensor, sep = ''))
    out <- list ('values' = eigen.dec $ values [1:n.tensor],
                 'matrices' = eigen.matrices)
    if (return.projection)
    {
      project <- aaply (matrix.array, 3, 
                        function (A, B) aaply (B, 3, frobenius.prob, y = A),
        B = eigen.matrices)
      out $ projection <- project
      out $ projection <- scale (out $ projection, center = FALSE, scale = TRUE)
    }
    return (out)
  }
