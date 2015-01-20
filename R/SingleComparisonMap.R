#' Generic Single Comparison Map functions for creating parallel list methods
#' Internal functions for making eficient comparisons.
#' @param matrix.list list of matrices being compared
#' @param y.mat single matrix being compared to list
#' @param MatrixCompFunc Function used to compare pair of matrices, must output a vector: comparisons and probabilities
#' @param ... Aditional arguments to MatrixCompFunc
#' @param num.cores Number of threads to use in computation. Requires doMC library.
#' @return Matrix of comparisons, matrix of probabilities.
#' @author Diogo Melo
#' @seealso \code{\link{MantelCor}}, \code{\link{KrzCor}},\code{\link{RandomSkewers}}
SingleComparisonMap  <- function(matrix.list, y.mat, MatrixCompFunc, ..., num.cores){
  if (num.cores > 1) {
    doMC::registerDoMC(num.cores)
    parallel = TRUE
  } else{
    parallel = FALSE
  }
  if(!all(c(laply(matrix.list, dim)) == dim(y.mat)[1]))
    stop("Matrices on list and single matrice dimension do not match")
  else
    out <- ldply(matrix.list,
                 function(x) {MatrixCompFunc(x, y.mat, ...)},
                 .parallel = parallel)
  return(out)
}
