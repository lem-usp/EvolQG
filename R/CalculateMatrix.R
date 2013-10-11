#' Calculate Covariance Matrix from LinearModel
#'
#' Calculates covariance matrix using the maximum likelihood estimator and the model residuals.
#' @param linear.m Linear model adjusted for original data.
#'
#' @value Estimated covariance matrix.
#' @references https://github.com/lem-usp/Morphometrics/wiki/
#' @author Diogo Melo, Fabio Machado
#'
#' @examples
#' data(iris)
#' options(contrasts=c("contr.sum","contr.poly"))
#' iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
#' cov.matrix <- CalculateMatrix(iris.lm)
#'
#' #To obtain a corrlation matrix, use:
#' cor.matrix <- cov2cor(cov.matrix)
#' @keyword covariancematrix
CalculateMatrix <- function(linear.m){
  cov.matrix = var(linear.m$residuals)*((dim(linear.m$residuals)[1]-1)/linear.m$df.residual)
  return (cov.matrix)
}
