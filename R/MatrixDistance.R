#' Matrix Riemann Distance
#'
#' Return distance between two covariance matrices
#'
#' @param cov.x covariance or correlation matrix
#' @param cov.y covariance or correlation matrix
#' @value Riemann distance between cov.x and cov.y
#' @export
#' @keywords matrixdistance
#' @keywords matrixcomparison

MatrixRiemannDist <- function(cov.x, cov.y) {
  return (sqrt(sum(log(eigen(solve(cov.x, cov.y))$values)^2)))
}
