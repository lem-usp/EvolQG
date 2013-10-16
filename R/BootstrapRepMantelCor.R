#' Bootstrap analysis via ressampling
#'
#'   Calculates the repeatability of the correlation matrix of the suplied data
#'   via bootstrap ressampling
#'
#'   Samples with replacement are taken from the full population, the correlation matrix is calculated
#'   and compared to the full population matrix via Mantel correlation 
#'
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param iterations Number of ressamples to take
#' @param num.cores Number of threads to use in computation. Requires doMC library.
#' @return returns the mean repeatability, or mean value of comparisons from samples to original statistic.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MonteCarloStat}}, \code{\link{AlphaRep}}, \code{\link{BootstrapRepKrzCor}}, \code{\link{BootstrapRep}}, \code{\link{BootstrapRepRandomSkewers}}
#' @export
#' @examples
#'
#' BootstrapRepMantelCor(iris[,1:4], 50)
#'
#'
#' @keywords bootstrap
#' @keywords repetabilities

BootstrapRepMantelCor <- function(ind.data, iterations = 1000, num.cores = 1){
  repeatability <- BootstrapRep(ind.data, iterations,
                                ComparisonFunc = function(x, y) MantelCor(x, y, 1)[1],
                                StatFunc = cor,
                                num.cores = num.cores)
  return(mean(repeatability))
}
