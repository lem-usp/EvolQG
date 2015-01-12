#' Bootstrap analysis via ressampling
#'
#'
#'   Calculates the repeatability of the covariance matrix of the suplied data
#'   via bootstrap ressampling
#'
#'   Samples with replacement are taken from the full population, a statistic calculated
#'   and compared to the full population statistic. 
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param iterations Number of ressamples to take
#' @param ComparisonFunc Comparison function for calculated statistic, either "randomskewers", "mantel" or "krzanowski" correlations
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. MantelCor always uses correlation matrix.
#' @param num.cores Number of threads to use in computation. Requires doMC library.
#' @return returns the mean repeatability, or mean value of comparisons from samples to original statistic.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MonteCarloStat}}, \code{\link{AlphaRep}}
#' @export
#' @examples
#' BootstrapRep(iris[,1:4], "mantel",
#'              iterations = 5,
#'              num.cores = 1)
#'              
#' BootstrapRep(iris[,1:4], "randomskewers", 50)
#'
#' # Partial matching for comparison function also works.
#' BootstrapRep(iris[,1:4], "krz", 50, TRUE)
#' 
#' @keywords bootstrap
#' @keywords repetabilities

BootstrapRep <- function(ind.data,
                         ComparisonFunc = c("randomskewers", "mantel", "krzanowski"),
                         iterations = 1000, 
                         correlation = F, 
                         num.cores = 1){
  ComparisonFunc = match.arg(ComparisonFunc)
  switch(ComparisonFunc,
        randomskewers = BootstrapRepRandomSkewers(ind.data, iterations, correlation, num.cores),
        mantel = BootstrapRepMantelCor(ind.data, iterations, num.cores),
        krzanowski = BootstrapRepKrzCor(ind.data, iterations, correlation, num.cores))
}

BootstrapRep_primitive <- function (ind.data, iterations,
                                    ComparisonFunc, StatFunc,
                                    num.cores){
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
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

BootstrapRepKrzCor <- function(ind.data, iterations = 1000, correlation = F, num.cores = 1){
  if(correlation)  StatFunc <- cor
  else             StatFunc <- cov
  repeatability <- BootstrapRep_primitive(ind.data, iterations,
                                          ComparisonFunc = KrzCor,
                                          StatFunc = StatFunc,
                                          num.cores = num.cores)
  return(mean(repeatability))
}

BootstrapRepMantelCor <- function(ind.data, iterations = 1000, num.cores = 1){
  repeatability <- BootstrapRep_primitive(ind.data, iterations,
                                          ComparisonFunc = function(x, y) MantelCor(x, y, 1)[1],
                                          StatFunc = cor,
                                          num.cores = num.cores)
  return(mean(repeatability))
}

BootstrapRepRandomSkewers <- function(ind.data, iterations = 1000, correlation = F, num.cores = 1){
  if(correlation)  StatFunc <- cor
  else             StatFunc <- cov
  repeatability <- BootstrapRep_primitive(ind.data, iterations,
                                          ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
                                          StatFunc = StatFunc,
                                          num.cores = num.cores)
  return(mean(repeatability))
}
