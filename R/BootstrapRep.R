#' Bootstrap analysis via resampling
#'
#'
#'   Calculates the repeatability of the covariance matrix of the suplied data
#'   via bootstrap resampling
#'
#'   Samples with replacement are taken from the full population, a statistic calculated
#'   and compared to the full population statistic. 
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param ComparisonFunc Comparison function for calculated statistic, either "randomskewers", "mantel" , "krzanowski" or "pcasimilarity" correlations
#' @param iterations Number of resamples to take
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. Mantel always uses correlation matrix.
#' @param num.cores Number of threads to use in computation.
#' The doMC library must be loaded.
#' @return returns the mean repeatability, that is, the mean value of comparisons from samples to original statistic.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MonteCarloStat}}, \code{\link{AlphaRep}}
#' @export
#' @examples
#' BootstrapRep(iris[,1:4], "mantel",
#'              iterations = 5)
#'              
#' BootstrapRep(iris[,1:4], "randomskewers", 50)
#'
#' # Partial matching for comparison function also works.
#' BootstrapRep(iris[,1:4], "krz", 50, TRUE)
#' 
#' #Multiple threads can be used with doMC library
#' library(doMC)
#' BootstrapRep(iris[,1:4], "mantel",
#'              iterations = 5,
#'              num.cores = 2)
#' @keywords bootstrap
#' @keywords repetabilities

BootstrapRep <- function(ind.data,
                         ComparisonFunc = c("randomskewers", 
                                            "mantel", 
                                            "krzanowski", 
                                            "pcasimilarity"),
                         iterations = 1000, 
                         correlation = FALSE, 
                         num.cores = 1){
  ComparisonFunc = match.arg(ComparisonFunc)
  switch(ComparisonFunc,
        randomskewers = BootstrapRepRandomSkewers(ind.data, iterations, correlation, num.cores),
        mantel = BootstrapRepMantelCor(ind.data, iterations, num.cores),
        krzanowski = BootstrapRepKrzCor(ind.data, iterations, correlation, num.cores),
        pcasimilarity = BootstrapRepPCAsimilarity(ind.data, iterations, correlation, num.cores))
}

BootstrapRep_primitive <- function (ind.data, iterations,
                                    ComparisonFunc, StatFunc,
                                    num.cores){
  if (num.cores > 1) {
    doMC::registerDoMC(num.cores)
    parallel = TRUE
  } else{
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

BootstrapRepKrzCor <- function(ind.data, iterations, correlation, num.cores){
  if(correlation)  StatFunc <- cor
  else             StatFunc <- cov
  repeatability <- BootstrapRep_primitive(ind.data, iterations,
                                          ComparisonFunc = KrzCor,
                                          StatFunc = StatFunc,
                                          num.cores = num.cores)
  return(mean(repeatability))
}

BootstrapRepMantelCor <- function(ind.data, iterations, num.cores){
  repeatability <- BootstrapRep_primitive(ind.data, iterations,
                                          ComparisonFunc = function(x, y) cor(x[lower.tri(x)], 
                                                                              y[lower.tri(y)]),
                                          StatFunc = cor,
                                          num.cores = num.cores)
  return(mean(repeatability))
}

BootstrapRepRandomSkewers <- function(ind.data, iterations, correlation, num.cores){
  if(correlation)  StatFunc <- cor
  else             StatFunc <- cov
  repeatability <- BootstrapRep_primitive(ind.data, iterations,
                                          ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
                                          StatFunc = StatFunc,
                                          num.cores = num.cores)
  return(mean(repeatability))
}

BootstrapRepPCAsimilarity <- function(ind.data, iterations, correlation, num.cores){
  if(correlation)  StatFunc <- cor
  else             StatFunc <- cov
  repeatability <- BootstrapRep_primitive(ind.data, iterations,
                                          ComparisonFunc = PCAsimilarity,
                                          StatFunc = StatFunc,
                                          num.cores = num.cores)
  return(mean(repeatability))
}

#' R2 confidence intervals by bootstrap resampling
#'
#' Random populations are generated by  resampling 
#' the suplied data or residuals. R2 is calculated on all the
#' random population's correlation matrices, provinding a distribution based on the original data.
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param iterations Number of resamples to take
#' @param num.cores Number of threads to use in computation.
#' The doMC library must be loaded.
#' @return returns a vector with the R2 for all populations
#' @author Diogo Melo Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
#' @export
#' @import plyr
#' @importFrom mvtnorm rmvnorm
#' @examples
#' r2.dist <- BootstrapR2(iris[,1:4], 30)
#' quantile(r2.dist)
#' @keywords bootstrap
#' @keywords integration
#' @keywords repeatability
BootstrapR2 <- function (ind.data, iterations = 1000, num.cores = 1) {
  it.r2 <- BootstrapRep_primitive(ind.data, iterations,
                          ComparisonFunc = function(x, y) y,
                          StatFunc = function(x) CalcR2(cor(x)),
                          num.cores = num.cores)
  return (it.r2)
}
