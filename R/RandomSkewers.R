#' Compare matrices via RandomSkewers
#'
#' Calculates covariance matrix correlation via random skewers
#'
#' @param cov.x Single covariance matrix or list of covariance matrices.
#' If single matrix is suplied, it is compared to cov.y.
#' If list is suplied and no cov.y is suplied, all matrices
#' are compared.
#' If cov.y is suplied, all matrices in list are compared to it.
#' @param cov.y First argument is compared to cov.y.
#' Optional if cov.x is a list.
#' @param iterations Number of random vectors used in comparison.
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param num.cores If list is passed, number of threads to use in computation. Requires doMC library.
#' @param ... aditional arguments passed to other methods
#' @return
#' If cov.x and cov.y are passed, returns average value
#' of random skewers ('AC'), significance ('prob') and standard deviation
#' of random skewers ('SD')
#'
#' If cov.x and cov.y are passed, same as above, but for all matrices in cov.x.
#'
#' If only a list is passed to cov.x, a matrix of RandomSkewers average
#' values and probabilities of all comparisons.
#' If repeat.vector is passed, comparison matrix is corrected below
#' diagonal and repeatabilities returned in diagonal.
#' @export
#' @rdname RandomSkewers
#' @references Cheverud, J. M., and Marroig, G. (2007). Comparing covariance matrices:
#' Random skewers method compared to the common principal components model.
#' Genetics and Molecular Biology, 30, 461-469.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{KrzCor}},\code{\link{MantelCor}}
#' @examples
#' c1 <- RandomMatrix(10)
#' c2 <- RandomMatrix(10)
#' c3 <- RandomMatrix(10)
#' RandomSkewers(c1, c2)
#'
#' RandomSkewers(list(c1, c2, c3))
#'
#' reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRepRandomSkewers, 10, 10))
#' RandomSkewers(list(c1, c2, c3), repeat.vector = reps)
#'
#' c4 <- RandomMatrix(10)
#' RandomSkewers(list(c1, c2, c3), c4)
#' @keyword matrixcomparison
#' @keyword matrixcorrelation
#' @keyword randomskewers
RandomSkewers <- function(cov.x, cov.y, ...) UseMethod("RandomSkewers")

#' @rdname RandomSkewers
#' @method RandomSkewers default
#' @S3method RandomSkewers default
RandomSkewers.default <- function (cov.x, cov.y, iterations = 1000, ...) {
  traits <- dim (cov.x) [1]
  base.vector <- Normalize(rnorm(traits))
  random.vectors <- array (rnorm (iterations * traits, mean = 0, sd = 1), c(traits, iterations))
  random.vectors <- apply (random.vectors, 2, Normalize)
  dist <- base.vector %*% random.vectors
  dz1 <- apply (cov.x %*% random.vectors, 2, Normalize)
  dz2 <- apply (cov.y %*% random.vectors, 2, Normalize)
  real <- apply (dz1 * dz2, 2, sum)
  ac <- mean (real)
  stdev <- sd (real)
  prob <- sum (ac < dist) / iterations
  output <- c(ac, prob, stdev)
  names(output) <- c("AC","P","SD")
  return(output)
}

#' @rdname RandomSkewers
#' @method RandomSkewers list
#' @S3method RandomSkewers list
RandomSkewers.list <- function (cov.x, cov.y = NULL, iterations = 1000, repeat.vector = NULL, num.cores = 1, ...)
{
  if (is.null (cov.y)) {
    out <- ComparisonMap(cov.x,
                         function(x, y) RandomSkewers.default(x, y, iterations),
                         repeat.vector = repeat.vector,
                         num.cores = num.cores)
  }
  else{
    out <- SingleComparisonMap(cov.x, cov.y,
                               function(x, y) RandomSkewers.default(x, y, iterations),
                               num.cores = num.cores)
  }
  return(out)
}
