#' Compare matrices via the correlation between response vectors
#'
#' Compares the expected response to selection for two matrices for a specific set of
#' selection gradients (not random gradients like in the RandomSkewers method)
#'
#' @param cov.x Single covariance matrix or list of covariance matrices.
#' If single matrix is suplied, it is compared to cov.y.
#' If list is suplied and no cov.y is suplied, all matrices
#' are compared.
#' If cov.y is suplied, all matrices in list are compared to it.
#' @param cov.y First argument is compared to cov.y.
#' Optional if cov.x is a list.
#' @param skewers matrix of column vectors to be used as gradients
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param ... aditional arguments passed to other methods.
#' @return vector of vector correlations between the expected responses for the two matrices for each supplied vector
#' @export
#' @useDynLib evolqg
#' @rdname DeltaZCorr
#' @references Cheverud, J. M., and Marroig, G. (2007). Comparing covariance matrices:
#' Random skewers method compared to the common principal components model.
#' Genetics and Molecular Biology, 30, 461-469.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{KrzCor}},\code{\link{MantelCor}}
#' @examples
#' x <- RandomMatrix(10, 1, 1, 10)
#' y <- RandomMatrix(10, 1, 1, 10)
#' 
#' n_skewers = 10
#' skewers = matrix(rnorm(10*n_skewers), 10, n_skewers)
#' DeltaZCorr(x, y, skewers)
#' 
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords randomskewers
DeltaZCorr <- function(cov.x, cov.y, skewers, ...) UseMethod("DeltaZCorr")

#' @rdname DeltaZCorr
#' @export
DeltaZCorr.default <- function (cov.x, cov.y, skewers, ...) {
  output <- as.numeric(delta_z_corr(cov.x, cov.y, ncol(skewers), skewers))
  names(output) = colnames(skewers)
  return(output)
}

#' @rdname DeltaZCorr
#' @method DeltaZCorr list
#' @export
DeltaZCorr.list <- function (cov.x, cov.y = NULL, skewers, parallel = FALSE, ...)
  {
    if (is.null (cov.y)) {
      output <- ComparisonMap(cov.x,
                              function(x, y) c(DeltaZCorr(x, y, skewers), NA),
                              parallel = parallel)
    } else{
      output <- SingleComparisonMap(cov.x, cov.y,
                                    function(x, y) c(DeltaZCorr(x, y, skewers), NA),
                                    parallel = parallel)
    }
    return(output)
  }