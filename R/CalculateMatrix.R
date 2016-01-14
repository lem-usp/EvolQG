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
  N = linear.m$df.residual
  p = dim(x)[2]
  if(is.null(nu)) nu = dim(x)[2]
  if(is.null(S_0)){
    S_0 = diag(0, p)
    diag(S_0) = diag(S_x/N) * nu
  }
  S_N <- S_0 + S_x
  nu_N <- nu + N
  MAP <- (S_0 + S_x) / (nu + N)
  MLE <- S_x/N
  if(!is.null(samples)){
   S_sample <- laply(rlply(samples, riwish(nu_N, S_N)), identity)
   median.P <- aaply(S_sample, 2:3, median)
   class(S_sample) <- "mcmc_sample"
   return(list(MAP = MAP, 
               MLE = MLE, 
               P = median.P, 
               Ps = S_sample))
  }
  else return(list(MAP = MAP, MLE = MLE))
}

#' @export
#' @method print mcmc_sample
#' @aliases CalculateMatrix
print.mcmc_sample <- function(x, ...) {
  cat(paste0("MCMC sample of ", dim(x)[1], 
             " matrices of dimension ", dim(x)[2], 
             "x", dim(x)[3]), "\n")
  cat("Median matrix:\n")
  print(aaply(x, 2:3, median))
}