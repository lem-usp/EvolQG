#' Calculate Mahalonabis distance for many vectors
#'
#' Calculates the Mahalanobis distance between a list of species mean, using a global covariance matrix
#' @param mean.list list of species means being compared
#' @param cov.matrix covariance matrix defining the metric tensor to be used
#' @param num.cores Number of threads to use in computation. The doMC library must be loaded.
#' @return returns a matrix of species-species distances.
#' @author Diogo Melo
#' @export
#' @seealso \code{\link{mahalanobis}}
#' @examples
#' mean.1 <- colMeans(matrix(rnorm(30*10), 30, 10))
#' mean.2 <- colMeans(matrix(rnorm(30*10), 30, 10))
#' mean.3 <- colMeans(matrix(rnorm(30*10), 30, 10))
#' mean.list <- list(mean.1, mean.2, mean.3)
#'
#' # If cov.matrix is identity, calculated distance is euclidian
#' euclidian <- MultiMahalanobis(mean.list, diag(rep(1, 10)))
#' # else, it is not
#' half.euclidian <- MultiMahalanobis(mean.list, diag(rep(0.5, 10)))
#' 
#' #Multiple threads can be used with doMC library
#' library(doMC)
#' MultiMahalanobis(mean.list, RandomMatrix(10), num.cores = 2)
MultiMahalanobis <- function (mean.list, cov.matrix, num.cores = 1) {
  if (num.cores > 1) {
    doMC::registerDoMC(num.cores)
    parallel = TRUE
  } else{
    parallel = FALSE
  }
  num.mean<- length(mean.list)
  if(is.null(names(mean.list))) {names(mean.list) <- 1:num.mean}
  mean.names <- names (mean.list)
  CompareToN <- function(n) ldply(mean.list[(n+1):num.mean],
                                  function(x) {mahalanobis(x, mean.list[[n]], cov.matrix)},
                                  .parallel = parallel)
  comparisons <- adply(1:(num.mean-1), 1,  CompareToN, .parallel = parallel)
  dists <- acast(comparisons[-4], X1~.id, value.var = 'V1')[,mean.names[-1]]
  distances <- array (0, c(num.mean, num.mean))
  distances[upper.tri(distances)] <- dists[upper.tri(dists, diag=T)]
  distances[lower.tri(distances)] <- t(distances)[lower.tri(distances)]
  rownames (distances) <- mean.names
  colnames (distances) <- mean.names
  return (distances)
}
