gamma <- function(S, S_0) sum(diag((x <- (S - S_0)) %*% t(x)))

#' Calculates Gamma-distance between modular hypothesis and covariance matrix
#' 
#' Combines and compares many modularity hypothesis to a covariance matrix. Comparison values are
#' ajusted to the number os zeros in the hypothesis using a linear regression. 
#' @param cor.matrix Correlation matrix
#' @param modularity.hypot Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @note Hypothesis can be named as column names, and these will be used to make labels in the 
#' output
#' @export
#' @references Marquez, E.J. 2008. A statistical framework for testing modularity in multidimensional data. Evolution 62:2688-2708. 
#' @examples
#' # Create a single modularity hypothesis:
#' hypot = rep(c(1, 0), each = 6)
#' cor.hypot = CreateHypotMatrix(hypot)
#' 
#' # First with an unstructured matrix:
#' un.cor = RandomMatrix(12)
#' MantelModTest(cor.hypot, un.cor)
#' 
#' # Now with a modular matrix:
#' modules = matrix(c(rep(c(1, 0, 0), each = 5),
#' rep(c(0, 1, 0), each = 5),
#' rep(c(0, 0, 1), each = 5)), 15)
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
ModularityGamma <- function(cor.matrix, modularity.hypot){
  cor.hypot.list = CombineHypot(modularity.hypot)
  if(is.list(cor.hypot.list)){
    if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- 1:dim(modularity.hypot)[2]
    raw.gamma = laply(cor.hypot.list, function(x) gamma(cor.matrix, cor.matrix * x))
    num.zeros = laply(cor.hypot.list, function(x) sum(x[lower.tri(x)] == 0))  
    gamma_df = data.frame(.id = names(cor.hypot.list), 
                          corrected.gamma = residuals(lm(raw.gamma ~ num.zeros)) + mean(raw.gamma),
                          raw.gamma, stringsAsFactors = FALSE)
    gamma_df = gamma_df[order(gamma_df$corrected.gamma),]
    rownames(gamma_df) <- NULL
    return(list(gamma_rank = gamma_df, "Modularity_hypothesis" = modularity.hypot[,strsplit(gamma_df[1,1], "_")[[1]]]))
  } 
  else gamma(cor.matrix, cor.matrix * cor.hypot.list)
}