#' Compare matrices using PCA similarity factor
#'
#' @param cov.x Single covariance matrix ou list of covariance matrices.
#' If cov.x is a single matrix, it is compared to cov.y.
#' If cov.x is a list and no cov.y is suplied, all matrices
#' are compared to each other.
#' If cov.x is a list and cov.y is suplied, all matrices in cov.x are compared to cov.y.
#' @param ... aditional arguments passed to other methods
#' @param cov.y First argument is compared to cov.y.
#' @param ret.dim number of retained dimensions in the comparison. Defaults to all.
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param num.cores If list is passed in cov.x, number of threads to use in computation.
#' The doMC library must be loaded.
#' @return Ratio of projected variance to total variance
#' @references Singhal, A. and Seborg, D. E. (2005), Clustering multivariate time-series data. J. Chemometrics, 19: 427-438. doi: 10.1002/cem.945
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
PCAsimilarity.default <- function(cov.x, cov.y, ret.dim = NULL, ...) {
  if (is.null(ret.dim))
    ret.dim = dim(cov.x)[1]

  eg.x <- eigen(cov.x)
  eg.y <- eigen(cov.y)
  eg.x.values <- eg.x$values[1:ret.dim]
  eg.y.values <- eg.y$values[1:ret.dim]
  eg.x.vectors <- eg.x$vectors[,1:ret.dim]
  eg.y.vectors <- eg.y$vectors[,1:ret.dim]

  total_var <- sum(eg.x.values %*% eg.y.values)

  return (c(PCAsimilarity = sum((eg.x.values %o% eg.y.values) * ((t(eg.x.vectors) %*% (eg.y.vectors))**2))/total_var))
}

#' @rdname PCAsimilarity
#' @method PCAsimilarity list
#' @export
PCAsimilarity.list <- function (cov.x, cov.y = NULL, ...,
                         repeat.vector = NULL, num.cores = 1) {
  if (is.null (cov.y)) {
    out <- ComparisonMap(cov.x,
                         function(x, y) return(c(PCAsimilarity(x, y, ...), NA)),
                         repeat.vector = repeat.vector,
                         num.cores = num.cores)
    out <- out[[1]]
  } else{
    out <- SingleComparisonMap(cov.x, cov.y,
                         function(x, y) return(c(PCAsimilarity(x, y, ...), NA)),
                               num.cores = num.cores)
    out <- out[,-length(out)]
  }
  return(out)
}
