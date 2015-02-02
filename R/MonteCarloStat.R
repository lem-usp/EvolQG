#' Parametric population samples with covariance or correlation matrices
#'
#' Using a multivariate normal model, random populations are generated
#' using the suplied covariance matrix. A statistic is calculated on the
#' random population and compared to the statistic calculated on the
#' original matrix. 
#'
#' @param cov.matrix Covariance matrix.
#' @param sample.size Size of the random populations
#' @param iterations Number of random populations
#' @param ComparisonFunc Comparison functions for the calculated statistic
#' @param StatFunc Function for calculating the statistic
#' @param num.cores If list is passed, number of threads to use in computation.
#' The doMC library must be loaded.
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
#' MonteCarloStat(cov.matrix, sample.size = 30, iterations = 100,
#'                ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
#'                StatFunc = cov,
#'                num.cores = 1)
#' #Multiple threads can be used with doMC library
#' library(doMC)
#' MonteCarloStat(cov.matrix, sample.size = 30, iterations = 100,
#'                ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
#'                StatFunc = cov,
#'                num.cores = 2)
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
    doMC::registerDoMC(num.cores)
    parallel = TRUE
  } else{
    parallel = FALSE
  }
  if(!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if(sum(diag(cov.matrix)) == dim(cov.matrix)[1]) warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling.")
  populations  <- alply(1:iterations, 1,
                        function(x) rmvnorm (sample.size, sigma = cov.matrix, method = 'chol'),
                        .parallel=parallel)
  comparisons <- ldply (populations, function (x) ComparisonFunc (cov.matrix, StatFunc(x)),
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
#' @param sample.size Size of the random populations.
#' @param ComparisonFunc comparison function.
#' @param ... Aditional arguments passed to ComparisonFunc.
#' @param iterations Number of random populations.
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. MantelCor and MatrixCor should always uses correlation matrix.
#' @param num.cores If list is passed, number of threads to use in computation.
#' The doMC library must be loaded.
#' @details Since this function uses multivariate normal model to generate populations, only covariance matrices should be used, even when computing repeatabilities for covariances matrices.
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
#' MonteCarloRep(cov.matrix, sample.size = 30, RandomSkewers, iterations = 50)
#' MonteCarloRep(cov.matrix, sample.size = 30, RandomSkewers, num.vectors = 100, 
#'               iterations = 50, correlation = TRUE)
#' MonteCarloRep(cov.matrix, sample.size = 30, MantelCor, permutations = 1, correlation = TRUE)
#' MonteCarloRep(cov.matrix, sample.size = 30, KrzCor)
#' MonteCarloRep(cov.matrix, sample.size = 30, KrzCor, correlation = TRUE)
#'
#' #Multiple threads can be used with doMC library
#' library(doMC)
#' MonteCarloRep(cov.matrix, 30, RandomSkewers, iterations = 100, num.cores = 2)
#'
#' #Creating repeatability vector for a list of matrices
#' mat.list <- RandomMatrix(10, 3, 1, 10)
#' laply(mat.list, MonteCarloRep, 30, KrzCor, correlation = TRUE)
#'
#' @keywords parametricsampling
#' @keywords montecarlo
#' @keywords repeatability
MonteCarloRep <- function(cov.matrix,
                          sample.size,
                          ComparisonFunc,
                          ...,
                          iterations = 1000, 
                          correlation = FALSE, 
                          num.cores = 1){
  if(correlation)  {StatFunc <- cov; c2v <- cov2cor
  } else {StatFunc <- cov; c2v <- function(x) x}
  repeatability <- MonteCarloStat(cov.matrix, sample.size, iterations,
                                  ComparisonFunc = function(x, y) ComparisonFunc(c2v(x), 
                                                                                 c2v(y), ...),
                                  StatFunc = StatFunc,
                                  num.cores = num.cores)
  return(mean(repeatability[,2]))
}

#' R2 confidence intervals by parametric sampling
#'
#' Using a multivariate normal model, random populations are generated
#' using the suplied covariance matrix. R2 is calculated on all the
#' random population, provinding a distribution based on the original matrix.
#'
#' @param cov.matrix Covariance matrix.
#' @param sample.size Size of the random populations
#' @param iterations Number of random populations
#' @param num.cores Number of threads to use in computation.
#' The doMC library must be loaded.
#' @details Since this function uses multivariate normal model to generate populations, only covariance matrices should be used.
#' @return returns a vector with the R2 for all populations
#' @author Diogo Melo Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @export
#' @import plyr
#' @importFrom mvtnorm rmvnorm
#' @examples
#' r2.dist <- MonteCarloR2(RandomMatrix(10, 1, 1, 10), 30)
#' quantile(r2.dist)
#' @keywords parametricsampling
#' @keywords montecarlo
#' @keywords repeatability
MonteCarloR2 <- function (cov.matrix, sample.size, iterations = 1000, num.cores = 1) {
  it.r2 <- MonteCarloStat(cov.matrix, sample.size, iterations,
                          ComparisonFunc = function(x, y) y,
                          StatFunc = function(x) CalcR2(cor(x)),
                          num.cores = num.cores)
  return (it.r2[,2])
}

