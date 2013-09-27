MonteCarloRep <- function (cov.matrix, sample.size, iterations,
                           ComparisonFunc, StatFunc,
                           num.cores = 1)
  # Calculates cov.matrix repeatability using parametric sampling
  #
  # Args:
  #     cov.matrix: covariance or correlation matrix.
  #               if cov.matrix is a correlation matrix will use MantelCor,
  #               else, will use RandomSkewers
  #     sample.size: number of sample.sizeivuals on each sample
  #     iterations: number of samples
  #     Comparisonfunc: Arbitrary function to compare 2 matrices. Must return single numeric value.
  # Return:
  #     mean correlation of sample covariance matrices with original input cov.matrix
{
  library(mvtnorm)
  library(plyr)
  library(reshape2)
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
  return (mean(comparisons))
}

MonteCarloRepRandomSkewers <- function(cov.matrix, sample.size, iterations = 1000, num.cores = 1){
  repeatability <- MonteCarloRep(cov.matrix, sample.size, iterations,
                                 ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
                                 StatFunc = cov,
                                 num.cores = num.cores)
  return(repeatability)
}

MonteCarloRepMantelCor <- function(cov.matrix, sample.size, iterations = 1000, num.cores = 1){
  repeatability <- MonteCarloRep(cov.matrix, sample.size, iterations,
                                 ComparisonFunc = function(x, y) MantelCor(x, y, 1)[1],
                                 StatFunc = function(x) cov2cor(cov(x)),
                                 num.cores = num.cores)
  return(repeatability)
}

MonteCarloRepKrzCor <- function(cov.matrix, sample.size, correlation = F, iterations = 1000, num.cores = 1){
  if(correlation)  StatFunc <- function(x) cov2cor(cov(x))
  else StatFunc <- cov
  repeatability <- MonteCarloRep(cov.matrix, sample.size, iterations,
                                 ComparisonFunc = function(x, y) MantelCor(x, y, 1)[1],
                                 StatFunc = StatFunc,
                                 num.cores = num.cores)
  return(repeatability)
}
