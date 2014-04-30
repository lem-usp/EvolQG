#' Matrix Riemann Distance
#'
#' Return distance between two covariance matrices
#'
#' @param cov.x
#' @param cov.y
#' @value Riemann distance between cov.x and cov.y
#' @export
#' @keywords matrixdistance
#' @keywords matrixcomparison

MatrixRiemannDist <- function(cov.x, cov.y) {
  return (sqrt(sum(log(eigen(solve(cov.x, cov.y))$values)^2)))
}
