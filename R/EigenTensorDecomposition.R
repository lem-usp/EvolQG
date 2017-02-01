#' Eigentensor Decomposition 
#' 
#' This function performs eigentensor decomposition on a set of covariance matrices.
#'
#' @param matrices k x k x m array of m covariance matrices with k traits;
#' @param return.projection Should we project covariance matrices into estimated eigentensors? Defaults to TRUE
#' @param ... aditional arguments for methods
#' @return List with the following components:
#' @return mean mean covariance matrices used to center the sample (obtained from \code{\link{MeanMatrix}})
#' @return mean.sqrt square root of mean matrix (saved for use in other functions, 
#' such as \code{\link{ProjectMatrix}} and \code{\link{RevertMatrix}})
#' @return values vector of ordered eigenvalues associated with eigentensors;
#' @return matrices array of eigentensor in matrix form;
#' @return projection matrix of unstandardized projected covariance matrices over eigentensors.
#' 
#' @details The number of estimated eigentensors is the minimum between the number of data points (m) and 
#' the number of independent variables (k(k + 1)/2) minus one, in a similar manner to the usual 
#' principal component analysis.
#' 
#' @references Basser P. J., Pajevic S. 2007. Spectral decomposition of a 4th-order 
#' covariance tensor: Applications to diffusion tensor MRI. Signal Processing. 87:220-236.
#' @references Hine E., Chenoweth S. F., Rundle H. D., Blows M. W. 2009. Characterizing 
#' the evolution of genetic variance using genetic covariance tensors. Philosophical 
#' transactions of the Royal Society of London. Series B, Biological sciences. 364:1567-78.
#' 
#' @author Guilherme Garcia, Diogo Melo
#' @seealso \code{\link{ProjectMatrix}}, \code{\link{RevertMatrix}}
#' 
#' @examples 
#' \dontrun{
#' data(dentus)
#'
#' dentus.vcv <- daply (dentus, .(species), function(x) cov(x[,-5]))
#'
#' dentus.vcv <- aperm(dentus.vcv, c(2, 3, 1))
#'
#' dentus.etd <- EigenTensorDecomposition(dentus.vcv, TRUE)
#'
#' # Plot some results
#' par(mfrow = c(1, 2))
#' plot(dentus.etd $ values, pch = 20, type = 'b', ylab = 'Eigenvalue')
#' plot(dentus.etd $ projection [, 1:2], pch = 20, 
#'      xlab = 'Eigentensor 1', ylab = 'Eigentensor 2')
#' text(dentus.etd $ projection [, 1:2],
#'      labels = rownames (dentus.etd $ projection), pos = 2)
#' 
#' # we can also deal with posterior samples of covariance matrices using plyr
#' 
#' dentus.models <- dlply(dentus, .(species), 
#'                        lm, formula = cbind(humerus, ulna, femur, tibia) ~ 1)
#'
#' dentus.matrices <- llply(dentus.models, BayesianCalculateMatrix, samples = 100)
#'
#' dentus.post.vcv <- laply(dentus.matrices, function (L) L $ Ps)
#'
#' dentus.post.vcv <- aperm(dentus.post.vcv, c(3, 4, 1, 2))
#'
#' # this will perform one eigentensor decomposition for each set of posterior samples
#' dentus.post.etd <- alply(dentus.post.vcv, 4, EigenTensorDecomposition)
#' 
#' # which would allow us to observe the posterior 
#' # distribution of associated eigenvalues, for instance
#' dentus.post.eval <- laply (dentus.post.etd, function (L) L $ values)
#' 
#' boxplot(dentus.post.eval, xlab = 'Index', ylab = 'Value', 
#'         main = 'Posterior Eigenvalue Distribution')
#' }
#' @importFrom matrixcalc frobenius.prod
#' @importFrom expm logm sqrtm
#'
#' @rdname EigenTensorDecomposition 
#' @export
#'
EigenTensorDecomposition <- function (matrices, return.projection = TRUE, ...) UseMethod('EigenTensorDecomposition')

#' @rdname EigenTensorDecomposition
#' @method EigenTensorDecomposition list
#' @export
EigenTensorDecomposition.list <- function (matrices, return.projection = TRUE, ...)
  {
    list2array <- laply(matrices, function(L) L)
    list2array <- aperm (list2array, c (2, 3, 1))
    EigenTensorDecomposition(list2array, return.projection, ...)
  }

#' @rdname EigenTensorDecomposition
#' @method EigenTensorDecomposition default
#' @export
EigenTensorDecomposition.default <- function (matrices, return.projection = TRUE, ...)
  {
    mean.matrix <- MeanMatrix(matrices, ...)
    mean.sqrt <- sqrtm(mean.matrix)
    mean.is <- solve(mean.sqrt)
    lce.array <- aaply(matrices, 3, 
                       function (A) logm(mean.is %*% A %*% mean.is))  
    lce.array <- aperm(lce.array, c(2, 3, 1))
    n.traits <- dim(lce.array) [1]
    n.matrix <- dim(lce.array) [3]
    
    ### Build Sigma
    variances <- aaply (matrices, 3, diag)
    covariances <- aaply (matrices, 3, function (x) x [lower.tri (x)])
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
      list (rownames (matrices),
            colnames (matrices),
            paste ('ET', 1:n.tensor, sep = ''))
    
    out <- list ('mean' = mean.matrix,
                 'mean.sqrt' = mean.sqrt,
                 'values' = eigen.dec $ values [1:n.tensor],
                 'matrices' = eigen.matrices)
   
    if (return.projection)
      out $ projection <- aaply (lce.array, 3, 
                          function (A, B) aaply (B, 3, frobenius.prod, y = A),
                          B = eigen.matrices)
    return (out)
  }