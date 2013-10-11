#' Bootstrap analysis via ressampling
#'
#'
#'   Calculates the repeatability of the covariance matrix of the suplied data
#'   via bootstrap ressampling
#'
#'   Samples with replacement are taken from the full population, the covariance matrix is calculated
#'   and compared to the full population matrix using Random Skewers. Prepackaged functions
#'   for common comparison functions and statistics are suplied.
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param iterations Number of ressamples to take
#' @param num.cores Number of threads to use in computation. Requires doMC library.
#' @return returns the mean repeatability, or mean value of comparisons from samples to original statistic.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MonteCarloStat}}, \code{\link{AlphaRep}}, \code{\link{BootstrapRepKrzCor}}, \code{\link{BootstrapRepMantelCor}}, \code{\link{BootstrapRep}}
#' @export
#' @examples
#'
#' BootstrapRepRandomSkewers(iris[,1:4], 50)
#'
#' #Creating repeatability vector for a list of matrices
#'
#' mat.list <- RandomMatrix(10, 3, 1, 10)
#' unlist(lapply(mat.list, BootstrapRepRandomSkewers), 30)
#'
#' @keywords bootstrap
#' @keywords repetabilities


BootstrapRepRandomSkewers <- function(ind.data, iterations = 1000, num.cores = 1){
  repeatability <- BootstrapRep(ind.data, iterations,
                                ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
                                StatFunc = cov,
                                num.cores = num.cores)
  return(mean(repeatability))
}
