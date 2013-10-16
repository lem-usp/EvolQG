#' Bootstrap analysis via ressampling
#'
#'
#'   Calculates the repeatability of the covariance matrix of the suplied data
#'   via bootstrap ressampling
#'
#'   Samples with replacement are taken from the full population, a statistic calculated
#'   and compared to the full population statistic. Prepackaged functions
#'   for common comparison functions and statistics are suplied.
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param iterations Number of ressamples to take
#' @param ComparisonFunc Comparison function for calculated statistic
#' @param StatFunc Statistic to be calculated
#' @param num.cores Number of threads to use in computation. Requires doMC library.
#' @return returns the mean repeatability, or mean value of comparisons from samples to original statistic.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MonteCarloStat}}, \code{\link{AlphaRep}}, \code{\link{BootstrapRepKrzCor}}, \code{\link{BootstrapRepMantelCor}}, \code{\link{BootstrapRepRandomSkewers}}
#' @export
#' @examples
#' BootstrapRep(iris[,1:4], iterations = 5,
#'              ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
#'              StatFunc = cov,
#'              num.cores = 1)
#'
#' @keywords bootstrap
#' @keywords repetabilities

BootstrapRep <- function (ind.data, iterations,
                          ComparisonFunc, StatFunc,
                          num.cores){
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else{
    parallel = FALSE
  }
  if(isSymmetric(as.matrix(ind.data))) stop("input appears to be a matrix, use residuals.")
  sample.size <-  dim (ind.data) [1]
  c.matrix <- StatFunc(ind.data)
  populations  <- alply(1:iterations, 1,
                        function(x) ind.data[sample (1:sample.size, sample.size, TRUE),],
                        .parallel = parallel)
  comparisons <- laply (populations, function (x) ComparisonFunc (c.matrix, StatFunc(x)),
                        .parallel = parallel)
  return (comparisons)
}

