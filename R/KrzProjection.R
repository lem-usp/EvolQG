#' Compare matrices via Modified Krzanowski Correlation
#'
#' Calculates the modified Krzanowski correlation between matrices,
#' projecting the variance in each principal  components of the first
#' matrix in to the ret.dim.2 components of the second for a list of matrices
#'
#' @param cov.x Single covariance matrix ou list of covariance matrices.
#' If single matrix is suplied, it is compared to cov.y.
#' If list is suplied and no cov.y is suplied, all matrices
#' are compared.
#' If cov.y is suplied, all matrices in list are compared to it.
#' @param ... aditional arguments passed to other methods
#' @param cov.y First argument is compared to cov.y.
#' Ignored if cov.marix.1 is a list.
#' @param ret.dim.1 number of retained dimensions for first matrix in the comparison,
#' @param ret.dim.2 number of retained dimensions for second matrix in the comparison,
#' default for nxn matrix is n/2-1
#' @param num.cores If list is passed, number of threads to use in computation.
#' Requires doMC library.
#' @param full.results if FALSE returns only total variance,
#' if TRUE also per PC variance.
#' @return Ratio of projected variance to total variance, and ratio of projectecd total in each PC
#' @references Krzanowski, W. J. (1979). Between-Groups Comparison of Principal
#' Components. Journal of the American Statistical Association, 74(367),
#' 703. doi:10.2307/2286995
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{RandomSkewers}},\code{\link{MantelCor}}
#' @rdname KrzProjection
#' @export
#' @import plyr
#' @importFrom reshape2 melt
#' @examples
#' c1 <- RandomMatrix(10)
#' c2 <- RandomMatrix(10)
#' c3 <- RandomMatrix(10)
#' KrzProjection(c1, c2)
#'
#' KrzProjection(list(c1, c2, c3))
#' KrzProjection(list(c1, c2, c3), 5, 4)
#' KrzProjection(list(c1, c2, c3), 4, 5)
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords Krzanowski

KrzProjection <- function(cov.x, cov.y, ...) UseMethod("KrzProjection")


#' @rdname KrzProjection
#' @method KrzProjection default
#' @S3method KrzProjection default
KrzProjection.default <- function (cov.x, cov.y, ret.dim.1 = NULL, ret.dim.2 = NULL, ...) {
  num.traits <- dim(cov.x)[1]
  if (is.null(ret.dim.1))
    ret.dim.1 <- num.traits/2 - 1
  if (is.null(ret.dim.2))
    ret.dim.2 <- num.traits/2 - 1
  eigen.cov.x <- eigen(cov.x)
  eVal.1 <- eigen.cov.x$values
  eVec.1 <- eigen.cov.x$vectors
  eVar.1 <- t(aaply(1:num.traits, 1, function(n) eVec.1[,n]*sqrt(eVal.1[n])))
  eVec.2 <- eigen(cov.y)$vectors
  SumSq  <- function(x) sum(x^2)
  MapProjection  <- function(x) SumSq(t(aaply(1:ret.dim.2, 1, function(n) eVar.1[,x]%*%eVec.2[,n])))
  ProjectionNorms <- aaply(1:ret.dim.1, 1, MapProjection)
  output <- list(total.variation = sum(ProjectionNorms)/sum(eVal.1),
                 per.PC = ProjectionNorms/eVal.1[1:ret.dim.1])
  return (output)
}

#' @rdname KrzProjection
#' @method KrzProjection list
#' @S3method KrzProjection list
KrzProjection.list <- function(cov.x, cov.y = NULL,
                               ret.dim.1 = NULL, ret.dim.2 = NULL,
                               num.cores = 1, full.results = FALSE, ...){
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else
    parallel = FALSE
  if(full.results){
    CompareToNProj <- function(n) llply(cov.x,
                                        function(x) {KrzProjection(x,
                                                                   cov.x[[n]],
                                                                   ret.dim.1, ret.dim.2)},
                                        .parallel = parallel)
  }
  else{
    CompareToNProj <- function(n) llply(cov.x[names(cov.x) != n],
                                        function(x) {KrzProjection(x,
                                                                   cov.x[[n]],
                                                                   ret.dim.1, ret.dim.2)[1]},
                                        .parallel = parallel)
  }
  if(is.null(names(cov.x))) {names(cov.x) <- 1:length(cov.x)}
  comparisons.proj <- llply(names(cov.x),
                             CompareToNProj,
                             .parallel = parallel)
  if(full.results){
    names(comparisons.proj) = names(cov.x)
    return(comparisons.proj)
  }
  else{
    comparisons.proj <- melt(comparisons.proj)
    comparisons.proj[,4] = names(cov.x)[(comparisons.proj[,4])]
    comparisons.proj = comparisons.proj[,-2]
    comparisons.proj = acast(comparisons.proj, L2~L1)[names(cov.x), names(cov.x)]
    diag(comparisons.proj) = 0.
  }
  return(comparisons.proj)
}
