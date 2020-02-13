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

#' Calculate Covariance Matrix from a linear model fitted with lm() using different estimators
#'
#' Calculates covariance matrix using the maximum likelihood estimator, the maximum a posteriori (MAP)
#' estimator under a regularized Wishart prior, and if the sample is large enough can give samples from the 
#' posterior and the median posterior estimator.
#' @param linear.m Linear model adjusted for original data
#' @param samples number os samples to be generated from the posterior. Requires sample size to be at least as large as the number of dimensions
#' @param ... aditional arguments, currently ignored
#' @param nu degrees of freedom in prior distribution, defaults to the number of traits (this can be a too strong prior)
#' @param S_0 cross product matrix of the prior. Default is to use the observed variances and zero covariance
#'
#' @return Estimated covariance matrices and posterior samples
#' @author Diogo Melo, Fabio Machado
#' @references Murphy, K. P. (2012). Machine learning: a probabilistic perspective. MIT press.
#' @references Schafer, J., e Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statistical applications in genetics and molecular biology, 4(1).
#' @export
#' @examples
#' data(iris)
#' iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
#' matrices <- BayesianCalculateMatrix(iris.lm, nu = 0.1, samples = 100)
#'
#' @keywords covariancematrix
#' @importFrom MCMCpack riwish
#' @importFrom stats median
BayesianCalculateMatrix <- function(linear.m, samples = NULL, ..., nu = NULL, S_0 = NULL){
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