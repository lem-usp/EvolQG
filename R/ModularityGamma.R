#' Calculates Gamma-distance between modular hypothesis and covariance matrix
#' 
#' Combines and compares many modularity hypothesis to a covariance matrix. Comparison values are
#' ajusted to the number os zeros in the hypothesis using a linear regression. Best hypothesis can
#' be assed using a jackknife procedure.
#' 
#' @param c.matrix Correlation or covariance matrix
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param modularity.hypot Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @param n number of jackknife resamplings
#' @param leave.out number of individuals to be left out of each jackknife, default is 20\%
#' @param ... aditional arguments to be passed to raply for the jackknife
#' @note Hypothesis can be named as column names, and these will be used to make labels in the 
#' output. 
#' 
#' Jackknife will return the best hypothesis for each sample.
#' @export
#' @aliases JackKnifeModularityGamma
#' @references Marquez, E.J. 2008. A statistical framework for testing modularity in multidimensional data. Evolution 62:2688-2708. 
#' @examples
#' # Creating a modular matrix:
#' modules = matrix(c(rep(c(1, 0, 0), each = 5),
#'                  rep(c(0, 1, 0), each = 5),
#'                  rep(c(0, 0, 1), each = 5)), 15)
#'                  
#' cor.hypot = CreateHypotMatrix(modules)[[4]]
#' hypot.mask = matrix(as.logical(cor.hypot), 15, 15)
#' mod.cor = matrix(NA, 15, 15)
#' mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
#' mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
#' diag(mod.cor) = 1
#' mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
#' 
#' # True hypothesis and a bunch of random ones.
#' hypothetical.modules = cbind(modules, matrix(sample(c(1, 0), 4*15, replace=TRUE), 15, 4))
#' 
#' # if hypothesis columns are not named they are assigned numbers
#' colnames(hypothetical.modules) <- letters[1:7]
#'  
#' ModularityGamma(mod.cor, hypothetical.modules)
#' 
#' random_var = runif(15, 1, 10)
#' mod.data = mvtnorm::rmvnorm(100, sigma = sqrt(outer(random_var, random_var)) * mod.cor)
#' out_jack = JackKnifeModularityGamma(mod.data, hypothetical.modules)
#' table(out_jack)
ModularityGamma <- function(c.matrix, modularity.hypot){
  cor.hypot.list = CombineHypot(modularity.hypot)
  if(is.list(cor.hypot.list)){
    if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- 1:dim(modularity.hypot)[2]
    gamma_df <- SingleModGamma(cor.hypot.list, c.matrix)
    return(list(gamma_rank = gamma_df, 
                "Modularity_hypothesis" = modularity.hypot[,strsplit(gamma_df[1,1], "_")[[1]]]))
  } 
  else gamma(c.matrix, c.matrix * cor.hypot.list)
}

gamma <- function(S, S_0) sum(diag((x <- (S - S_0)) %*% t(x)))

SingleModGamma <- function(cor.hypot.list, c.matrix) {
  raw.gamma = laply(cor.hypot.list, function(x) gamma(c.matrix, c.matrix * x))
  num.zeros = laply(cor.hypot.list, function(x) sum(x[lower.tri(x)] == 0))  
  gamma_df = data.frame(.id = names(cor.hypot.list), 
                        corrected.gamma = residuals(lm(raw.gamma ~ num.zeros)) + mean(raw.gamma),
                        raw.gamma, stringsAsFactors = FALSE)
  gamma_df = gamma_df[order(gamma_df$corrected.gamma),]
  rownames(gamma_df) <- NULL
  return(gamma_df)
}

#' @export
#' @rdname ModularityGamma
JackKnifeModularityGamma <- function(ind.data, modularity.hypot, 
                                     n = 20, leave.out = floor(dim(ind.data)[1]/20),
                                     ...){
  if(isSymmetric(as.matrix(ind.data))) stop("input appears to be a matrix, use residuals")
  cor.hypot.list = CombineHypot(modularity.hypot)
  if(is.list(cor.hypot.list)){
    n.ind = dim(ind.data)[1]
    if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- 1:dim(modularity.hypot)[2]
    raply(n,function() {
      current.sample = sample(1:n.ind, n.ind - leave.out)
      SingleModGamma(cor.hypot.list, cov(ind.data[current.sample,]))[1,1]
    }, ...)
  }
  else stop("use at least two hypothesis vectors, preferably more")
}
