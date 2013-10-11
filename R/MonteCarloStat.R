#' Parametric repeatabilities with covariance or correlation matrices
#'
#' Using a multivariate normal model, random populations are generated
#' using the suplied covariance matrix. A statistic is calculated on the
#' random population and compared to the statistic calculated on the
#' original matrix. Prepackaged functions for common comparison functions
#' and statistics are suplied.
#'
#' @aliases MonteCarloRepRandomSkewers MonteCarloRepMantelCor MonteCarloRepKrzCor MonteCarloR2
#' @param cov.matrix Covariance matrix.
#' @param sample.size Size of the random populations
#' @param iterations Number of random populations
#' @param ComparisonFunc Comparison functions for the calculated statistic
#' @param StatFunc Function for calculating the statistic
#' @param correlation When using BootstrapRepKrzCor, statistic can be correlation or covariance. If TRUE, uses correlation.
#' @param num.cores Number of threads to use in computation. Requires doMC library.
#' @details Since this function uses multivariate normal model to generate populations, only covariance matrices should be used.
#' @return returns the mean repeatability, or mean value of comparisons from samples to original statistic.
#' @author Diogo Melo Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @rdname MonteCarloStat
#' @export
#' @import plyr
#' @importFrom mvtnorm rmvnorm
#' @examples
#' cov.matrix <- RandomMatrix(10, 1, 1, 10)
#'
#' MonteCarloStat(cov.matrix, sample.size = 30, iterations = 1000,
#'                ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
#'                StatFunc = cov,
#'                num.cores = 1)
#'
#'
#' MonteCarloRepRandomSkewers(cov.matrix, 30)
#' MonteCarloRepMantelCor(cov.matrix, 30)
#' MonteCarloRepKrzCor(cov.matrix, 30)
#' MonteCarloRepKrzCor(cov.matrix, 30, TRUE)
#'
#' #Creating repeatability vector for a list of matrices
#' mat.list <- RandomMatrix(10, 3, 1, 10)
#' unlist(lapply(mat.list, MonteCarloRepKrzCor, 30, TRUE))
#'
#' #Calculating R2 confidence intervals
#' r2.dist <- MonteCarloR2(RandomMatrix(10, 1, 1, 10), 30)
#' quantile(r2.dist)
#' @keywords parametricsampling
#' @keywords montecarlo
#' @keywords repeatability
MonteCarloStat <- function (cov.matrix, sample.size, iterations,
                            ComparisonFunc, StatFunc,
                            num.cores = 1) {
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else{
    parallel = FALSE
  }
  if(sum(diag(cov.matrix)) == dim(cov.matrix)[1]) warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling.")
  populations  <- alply(1:iterations, 1,
                        function(x) rmvnorm (sample.size, sigma = cov.matrix, method = 'chol'),
                        .parallel=parallel)
  comparisons <- laply (populations, function (x) ComparisonFunc (cov.matrix, StatFunc(x)),
                        .parallel = parallel)
  return (comparisons)
}

#' @export
#' @rdname MonteCarloStat
MonteCarloRepRandomSkewers <- function(cov.matrix, sample.size, iterations = 1000, num.cores = 1){
  repeatability <- MonteCarloStat(cov.matrix, sample.size, iterations,
                                  ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
                                  StatFunc = cov,
                                  num.cores = num.cores)
  return(mean(repeatability))
}

#' @export
#' @rdname MonteCarloStat
MonteCarloRepMantelCor <- function(cov.matrix, sample.size, iterations = 1000, num.cores = 1){
  repeatability <- MonteCarloStat(cov.matrix, sample.size, iterations,
                                  ComparisonFunc = function(x, y) MantelCor(x, y, 1)[1],
                                  StatFunc = function(x) cov2cor(cov(x)),
                                  num.cores = num.cores)
  return(mean(repeatability))
}

#' @export
#' @rdname MonteCarloStat
MonteCarloRepKrzCor <- function(cov.matrix, sample.size, correlation = F, iterations = 1000, num.cores = 1){
  if(correlation)  StatFunc <- function(x) cov2cor(cov(x))
  else StatFunc <- cov
  repeatability <- MonteCarloStat(cov.matrix, sample.size, iterations,
                                  ComparisonFunc = function(x, y) MantelCor(x, y, 1)[1],
                                  StatFunc = StatFunc,
                                  num.cores = num.cores)
  return(mean(repeatability))
}

#' @export
#' @rdname MonteCarloStat
MonteCarloR2 <- function (cov.matrix, sample.size, iterations = 1000, num.cores = 1) {
  it.r2 <- MonteCarloStat(cov.matrix, sample.size, iterations,
                          ComparisonFunc = function(x, y) y,
                          StatFunc = function(x) CalcR2(cor(x)),
                          num.cores = num.cores)
  return (it.r2)
}

