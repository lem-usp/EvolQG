#' Compare matrices via Mantel Correlation
#'
#' Calculates correlation matrix correlation via Mantel.
#'
#' @param cor.x Single covariance matrix or list of covariance matrices.
#'
#' If single matrix is suplied, it is compared to cor.y.
#'
#' If list is suplied and no cor.y is suplied, all matrices
#' are compared.
#'
#' If cor.y is suplied, all matrices in list are compared to it.
#' @param cor.y First argument is compared to cor.y.
#' Optional if cor.x is a list.
#' @param iterations Number of permutations used in significance calculation.
#' @param mod Set TRUE to use mantel in testing modularity hipotesis. Will return
#' AVG+, AVG- and AVG Ratio based on binary hipotesis matrix.
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param num.cores If list is passed, number of threads to use in computation. Requires doMC library.
#' @param ... aditional arguments passed to other methods
#' @return If cor.x and cor.y are passed, returns matrix pearson
#' correlation and significance via mantel permutations.
#'
#' If cor.x and cor.y are passed, same as above, but for all matrices in cor.x.
#'
#' If only cor.x is passed, a matrix of MantelCor average
#' values and probabilities of all comparisons.
#' If repeat.vector is passed, comparison matrix is corrected above
#' diagonal and repeatabilities returned in diagonal.
#' @export
#' @importFrom vegan mantel
#' @rdname MantelCor
#' @references http://en.wikipedia.org/wiki/Mantel_test
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{KrzCor}},\code{\link{RandomSkewers}},\code{\link{mantel}},\code{\link{RandomSkewers}},\code{\link{TestModularity}}
#' @examples
#' c1 <- RandomMatrix(10)
#' c2 <- RandomMatrix(10)
#' c3 <- RandomMatrix(10)
#' MantelCor(c1, c2)
#'
#' MantelCor(list(c1, c2, c3))
#'
#' reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRepMantelCor, 10, 10))
#' MantelCor(list(c1, c2, c3), repeat.vector = reps)
#'
#' c4 <- RandomMatrix(10)
#' MantelCor(list(c1, c2, c3), c4)
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords randomskewers
MantelCor <- function (cor.x, cor.y, ...) UseMethod("MantelCor")

#' @rdname MantelCor
#' @method MantelCor default
#' @S3method MantelCor default
MantelCor.default <- function (cor.x, cor.y, iterations = 1000, mod = FALSE, ...) {
  mantel.out <- mantel(cor.x, cor.y, permutations = iterations)
  correlation <- mantel.out$statistic
  prob <- mantel.out$signif
  if (mod == TRUE){
    index <- cor.y[lower.tri(cor.y)]
    avg.plus <- mean (cor.x [lower.tri(cor.x)] [index != 0])
    avg.minus <- mean (cor.x [lower.tri(cor.x)] [index == 0])
    avg.ratio <- avg.plus / avg.minus
    output <- c(correlation, prob, avg.plus, avg.minus, avg.ratio)
    names(output) <- c("Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio")
  }
  else{
    if(sum(diag(cor.x)) != dim(cor.x)[1] | sum(diag(cor.y))!= dim(cor.y)[1])
      warning("Matrices do not appear to be correlation matrices. Use with caution.")
    output <- c(correlation, prob)
    names(output) <- c("Rsquared", "Probability")
  }
  return (output)
}

#' @rdname MantelCor
#' @method MantelCor list
#' @S3method MantelCor list
MantelCor.list <- function (cor.x, cor.y = NULL,
                            iterations = 1000, repeat.vector = NULL,
                            mod = FALSE, num.cores = 1, ...)
{
  if (is.null (cor.y)) {
    out <- ComparisonMap(cor.x,
                         function(x, cor.y) MantelCor(x, cor.y, iterations),
                         repeat.vector = repeat.vector,
                         num.cores = num.cores)
  }
  else{
    out <- SingleComparisonMap(cor.x, cor.y,
                               function(x, y) MantelCor(y,
                                                                         x,
                                                                         iterations, mod = mod),
                               num.cores = num.cores)
  }
  return(out)
}
