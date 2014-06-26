#' Compare matrices using PCA similarity factor
#'
#' @param cov.x Single covariance matrix ou list of covariance matrices.
#' If single matrix is suplied, it is compared to cov.y.
#' If list is suplied and no cov.y is suplied, all matrices
#' are compared.
#' If cov.y is suplied, all matrices in list are compared to it.
#' @param ... aditional arguments passed to other methods
#' @param cov.y First argument is compared to cov.y.
#' If cov.x is a list, every element in cov.x is projected in cov.y.
#' @param ret.dim number of retained dimensions for matrix comparison,
#' default for nxn matrix is n/2-1
#' @param num.cores If list is passed, number of threads to use in computation.
#' Requires doMC library.
#' @return Ratio of projected variance to total variance
#' @references Krzanowski, W. J. (1979). Between-Groups Comparison of Principal
#' Components. Journal of the American Statistical Association, 74(367),
#' 703. doi:10.2307/2286995
#' @author Edgar Zanella Alvarenga
#' @seealso \code{\link{KrzProjection}},\code{\link{KrzCor}},\code{\link{RandomSkewers}},\code{\link{MantelCor}}
#' @rdname PCAsimilarity
#' @export
#' @import plyr
#' @examples
#' c1 <- RandomMatrix(10)
#' c2 <- RandomMatrix(10)
#' PCAsimilarity(c1, c2)
#'
#' m.list <- RandomMatrix(10, 3)
#' PCAsimilarity(m.list)
#' PCAsimilarity(m.list, ret.dim = 5)
#'
#' PCAsimilarity(m.list, c1)
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords Krzanowski
#' @keywords PCA

PCAsimilarity <- function(cov.x, cov.y, ...) UseMethod("PCAsimilarity")

#' @rdname PCAsimilarity
#' @method PCAsimilarity default
#' @export
PCAsimilarity.default <- function(cov.x, cov.y, ret.dim=NULL) {
  if (is.null(ret.dim)) {
    ret.dim = round(dim(cov.x)[1]/2 - 1)
  }
  eg.x <- eigen(cov.x)
  eg.y <- eigen(cov.y)
  total_var <- sum(eg.x$values %*% eg.y$values)

  return (sum((eg.x$values %o% eg.y$values) * ((t(eg.x$vectors) %*% (eg.y$vectors))**2))/total_var)
}

#' @rdname PCAsimilarity
#' @method PCAsimilarity list
#' @export
PCAsimilarity.list <- function (cov.x, cov.y = NULL,
                         ret.dim = NULL, repeat.vector = NULL,
                         num.cores = 1, ...) {
  if (is.null (cov.y)) {
    out <- ComparisonMap(cov.x,
                         function(x, y) return(c(PCAsimilarity(x, y, ret.dim), NA)),
                         repeat.vector = repeat.vector,
                         num.cores = num.cores)
    out <- out[[1]]
  } else{
    out <- SingleComparisonMap(cov.x, cov.y,
                         function(x, y) return(c(PCAsimilarity(x, y, ret.dim), NA)),
                               num.cores = num.cores)
    out <- out[,-length(out)]
  }
  return(out)
}
