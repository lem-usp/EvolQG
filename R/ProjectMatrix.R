#' Project Covariance Matrix
#'
#' This function projects a given covariance matrix into the basis provided by an eigentensor decomposition.
#'
#' @param matrix A symmetric covariance matrix for k traits
#' @param etd Eigentensor decomposition of m covariance matrices for k traits
#' (obtained from \code{\link{EigenTensorDecomposition}})
#' 
#' @return Vector of scores of given covariance matrix onto eigentensor basis.
#' 
#' @references Basser P. J., Pajevic S. 2007. Spectral decomposition of a 4th-order 
#' covariance tensor: Applications to diffusion tensor MRI. Signal Processing. 87:220-236.
#' @references Hine E., Chenoweth S. F., Rundle H. D., Blows M. W. 2009. Characterizing 
#' the evolution of genetic variance using genetic covariance tensors. Philosophical 
#' transactions of the Royal Society of London. Series B, Biological sciences. 364:1567-78.
#' 
#' @author Guilherme Garcia, Diogo Melo
#' @seealso \code{\link{EigenTensorDecomposition}}, \code{\link{RevertMatrix}}
#' 
#' @examples
#' # this function is useful for projecting posterior samples for a set of 
#' # covariance matrices onto the eigentensor decomposition done 
#' # on their estimated means
#' \dontrun{
#' data(dentus)
#' 
#' dentus.models <- dlply(dentus, .(species), lm, 
#'                        formula = cbind(humerus, ulna, femur, tibia) ~ 1)
#'
#' dentus.matrices <- llply(dentus.models, BayesianCalculateMatrix, samples = 100)
#'
#' dentus.post.vcv <- laply(dentus.matrices, function (L) L $ Ps)
#' dentus.post.vcv <- aperm(dentus.post.vcv, c(3, 4, 1, 2))
#' 
#' dentus.mean.vcv <- aaply(dentus.post.vcv, 3, MeanMatrix)
#' dentus.mean.vcv <- aperm(dentus.mean.vcv, c(2, 3, 1))
#' 
#' dentus.mean.etd <- EigenTensorDecomposition(dentus.mean.vcv)
#' dentus.mean.proj <- data.frame('species' = LETTERS [1:5], dentus.mean.etd $ projection)
#' 
#' dentus.post.proj <- adply(dentus.post.vcv, c(3, 4), ProjectMatrix, etd = dentus.mean.etd)
#' colnames(dentus.post.proj) [1:2] <- c('species', 'sample')
#' levels(dentus.post.proj $ species) <- LETTERS[1:5]
#'
#' require(ggplot2)
#' ggplot() +
#'   geom_point(aes(x = ET1, y = ET2, color = species), 
#'      data = dentus.mean.proj, shape = '+', size = 8) +
#'   geom_point(aes(x = ET1, y = ET2, color = species), 
#'      data = dentus.post.proj, shape = '+', size = 3) +
#'   theme_bw()
#' }
#' @importFrom matrixcalc frobenius.prod
#' @importFrom expm logm sqrtm
#'
#' @export
#' 
ProjectMatrix <-
  function (matrix, etd)
  {
    mean.is <- solve (etd $ mean.sqrt)
    log.center.mat <- logm (mean.is %*% matrix %*% mean.is)
    projection <- aaply (etd $ matrices, 3, frobenius.prod, y = log.center.mat)
    projection
  }
