#' Calculate Covariance Matrix from a linear model fitted with lm()
#'
#' Calculates covariance matrix using the maximum likelihood estimator and the model residuals.
#' @param linear.m Linear model adjusted for original data.
#'
#' @return Estimated covariance matrix.
#' @references https://github.com/lem-usp/evolqg/wiki/
#' @author Diogo Melo, Fabio Machado
#' @export
#' @examples
#' data(iris)
#' options(contrasts=c("contr.sum","contr.poly"))
#' iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
#' cov.matrix <- CalculateMatrix(iris.lm)
#'
#' #To obtain a corrlation matrix, use:
#' cor.matrix <- cov2cor(cov.matrix)
#' @keywords covariancematrix
CalculateMatrix <- function(linear.m){
  cov.matrix = var(linear.m$residuals)*((dim(linear.m$residuals)[1]-1)/linear.m$df.residual)
  return (cov.matrix)
}

#' @importFrom MCMCpack riwish
CalculateMatrix_Baysean <- function(linear.m, samples = NULL, ..., nu = NULL, S_0 = NULL){
  x <- linear.m$residuals
  S_x = t(x) %*% x
  N = dim(x)[1]
  p = dim(x)[2]
  if(is.null(nu)) nu = dim(x)[2] + 1
  if(is.null(S_0)){
    S_0 = diag(0, p)
    diag(S_0) = diag(S_x/N) * nu
  }
  S_N <- S_0 + S_x
  nu_N <- nu + N
  MAP <- (S_0 + S_x) / (nu + N)
  if(!is.null(samples)){
   S_sample <- rlply(samples, riwish(nu_N, S_N))
   return(list(MAP = MAP, posterior_samples = S_sample))
  }
  else return(MAP)
}