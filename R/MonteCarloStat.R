#' Parametric population samples with covariance or correlation matrices
#'
#' Using a multivariate normal model, random populations are generated
#' using the suplied covariance matrix. A statistic is calculated on the
#' random population and compared to the statistic calculated on the
#' original matrix. 
#'
#' @aliases MonteCarloR2
#' @param cov.matrix Covariance matrix.
#' @param sample.size Size of the random populations
#' @param iterations Number of random populations
#' @param ComparisonFunc Comparison functions for the calculated statistic
#' @param StatFunc Function for calculating the statistic
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

#' Parametric repeatabilities with covariance or correlation matrices
#'
#' Using a multivariate normal model, random populations are generated
#' using the suplied covariance matrix. A statistic is calculated on the
#' random population and compared to the statistic calculated on the
#' original matrix. 
#'
#' @param cov.matrix Covariance matrix.
#' @param ComparisonFunc Comparison function for calculated statistic, either "randomskewers", "mantel" or "krznowski" correlations
#' @param sample.size Size of the random populations
#' @param iterations Number of random populations
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. MantelCor always uses correlation matrix.
#' @param num.cores Number of threads to use in computation. Requires doMC library.
#' @details Since this function uses multivariate normal model to generate populations, only covariance matrices should be used.
#' @return returns the mean repeatability, or mean value of comparisons from samples to original statistic.
#' @author Diogo Melo Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @rdname MonteCarloRep
#' @export
#' @import plyr
#' @importFrom mvtnorm rmvnorm
#' @examples
#' cov.matrix <- RandomMatrix(10, 1, 1, 10)
#'
#' MonteCarloRep(cov.matrix, "randomskewers", 30)
#' MonteCarloRep(cov.matrix, "mantel", 30)
#' MonteCarloRep(cov.matrix, "krz", 30)
#' MonteCarloRep(cov.matrix, "krz", 30, TRUE)
#'
#' #Creating repeatability vector for a list of matrices
#' mat.list <- RandomMatrix(10, 3, 1, 10)
#' unlist(lapply(mat.list, MonteCarloRep, "krz", 30, correlation = TRUE))
#'
#' @keywords parametricsampling
#' @keywords montecarlo
#' @keywords repeatability
MonteCarloRep <- function(cov.matrix,
                          ComparisonFunc = c("randomskewers", "mantel", "krzanowski"),
                          sample.size,
                          iterations = 1000, 
                          correlation = F, 
                          num.cores = 1){
  ComparisonFunc = match.arg(ComparisonFunc)
  switch(ComparisonFunc,
         randomskewers = MonteCarloRepRandomSkewers(cov.matrix, sample.size, iterations, correlation, num.cores),
         mantel = MonteCarloRepMantelCor(cov.matrix, sample.size, iterations, num.cores),
         krzanowski = MonteCarloRepKrzCor(cov.matrix, sample.size, iterations, correlation, num.cores))
}

MonteCarloRepRandomSkewers <- function(cov.matrix, sample.size, iterations = 1000, correlation, num.cores = 1){
  if(correlation)  StatFunc <- cor
  else             StatFunc <- cov
  repeatability <- MonteCarloStat(cov.matrix, sample.size, iterations,
                                  ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
                                  StatFunc = cov,
                                  num.cores = num.cores)
  return(mean(repeatability))
}

MonteCarloRepMantelCor <- function(cov.matrix, sample.size, iterations = 1000, num.cores = 1){
  repeatability <- MonteCarloStat(cov.matrix, sample.size, iterations,
                                  ComparisonFunc = function(x, y) MantelCor(cov2cor(x), y, 1)[1],
                                  StatFunc = cor,
                                  num.cores = num.cores)
  return(mean(repeatability))
}

MonteCarloRepKrzCor <- function(cov.matrix, sample.size, correlation = F, iterations = 1000, num.cores = 1){
  if(correlation)  {StatFunc <- cor; c2v <- cov2cor}
  else {StatFunc <- cov; c2v <- function(x) x}
  repeatability <- MonteCarloStat(cov.matrix, sample.size, iterations,
                                  ComparisonFunc = function(x, y) KrzCor(c2v(x), c2v(y), 1)[1],
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

