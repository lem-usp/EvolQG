#' Matrix Riemann Distance
#'
#' Return distance between two covariance matrices
#'
#' @param cov.x covariance or correlation matrix
#' @param cov.y covariance or correlation matrix
#' @return Riemann distance between cov.x and cov.y
#' @export
#' @keywords matrixdistance
#' @keywords matrixcomparison

MatrixRiemannDist <- function(cov.x, cov.y) {
  return (sqrt(sum(log(eigen(solve(cov.x, cov.y))$values)^2)))
}

#'Distribution overlap distance
#'
#'Calculates the overlap between two normal distributions, 
#'defined as the probability that a draw from one distribution 
#'comes from the other
#'
#' @param cov.x covariance or correlation matrix
#' @param cov.y covariance or correlation matrix
#' @param iterations number of drows
#' @return Overlap distance between cov.x and cov.y
#'
#' @export
#' @importFrom mvtnorm rmvnorm dmvnorm
MatrixOverlapDist <- function(cov.x, cov.y, iterations = 10000){
  x_samples <- rmvnorm(iterations, sigma = cov.x)
  y_density <- dmvnorm(x_samples, sigma = cov.y)
  x_density <- dmvnorm(x_samples, sigma = cov.x)
  overlap <- mean(y_density/(x_density + y_density))  
  return(sqrt(1-overlap*2))
}
