#' Compare matrices via Krzanowski Correlation
#'
#' Calculates covariance matrix correlation via Krzanowski Correlation
#'
#' @param cov.x Single covariance matrix or list of covariance matrices.
#' If single matrix is suplied, it is compared to cov.y.
#' If list is suplied and no cov.y is suplied, all matrices
#' are compared.
#' If cov.y is suplied, all matrices in list are compared to it.
#' @param cov.y First argument is compared to cov.y.
#' Optional if cov.x is a list.
#' @param ret.dim number of retained dimensions in the comparison,
#' default for nxn matrix is n/2-1
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param num.cores If list is passed, number of threads to use in computation. Requires doMC library.
#' @param ... aditional arguments passed to other methods
#' @return If cov.x and cov.y are passed, returns Kzranowski correlation
#' 
#' If cov.x is a list and cov.y is passed, same as above, but for all matrices in cov.x.
#'
#' If only a list is passed to cov.x, a matrix of Kzranowski correlation 
#' values.
#' If repeat.vector is passed, comparison matrix is corrected above
#' diagonal and repeatabilities returned in diagonal.
#' @export
#' @rdname KrzCor
#' @references Krzanowski, W. J. (1979). Between-Groups Comparison of Principal
#' Components. Journal of the American Statistical Association, 74(367),
#' 703. doi:10.2307/2286995
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{RandomSkewers}},\code{\link{KrzProjection}},\code{\link{MantelCor}}
#' @examples
#' c1 <- RandomMatrix(10)
#' c2 <- RandomMatrix(10)
#' c3 <- RandomMatrix(10)
#' KrzCor(c1, c2)
#'
#' KrzCor(list(c1, c2, c3))
#'
#' reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRep, "krz", 10, 10))
#' KrzCor(list(c1, c2, c3), repeat.vector = reps)
#'
#' c4 <- RandomMatrix(10)
#' KrzCor(list(c1, c2, c3), c4)
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords Krzanowski
KrzCor <- function (cov.x, cov.y, ...) UseMethod("KrzCor")

#' @rdname KrzCor
#' @method KrzCor default
#' @export
KrzCor.default <- function (cov.x, cov.y, ret.dim = NULL, ...) {
  if (is.null(ret.dim))
    ret.dim = round(dim(cov.x)[1]/2 - 1)

  eg.x <- eigen(cov.x)
  eg.y <- eigen(cov.y)

  return (sum((t(eg.x$vectors[,1:ret.dim]) %*% (eg.y$vectors[,1:ret.dim]))**2)/ret.dim)
}

#' @rdname KrzCor
#' @method KrzCor list
#' @export
KrzCor.list <- function (cov.x, cov.y = NULL,
                         ret.dim = NULL, repeat.vector = NULL,
                         num.cores = 1, ...) {
  if (is.null (cov.y)) {
    out <- ComparisonMap(cov.x,
                         function(x, y) return(c(KrzCor(x, y, ret.dim), NA)),
                         repeat.vector = repeat.vector,
                         num.cores = num.cores)
    out <- out[[1]]
  } else{
    out <- SingleComparisonMap(cov.x, cov.y,
                         function(x, y) return(c(KrzCor(x, y, ret.dim), NA)),
                               num.cores = num.cores)
    out <- out[,-length(out)]
  }
  return(out)
}
