#' Modularity and integration analysis tool
#' 
#' Combines and compares many modularity hypothesis to a covariance matrix. Comparison values are
#' ajusted to the number os zeros in the hypothesis using a linear regression. Best hypothesis can
#' be assessed using a jack-knife procedure.
#' 
#' @param c.matrix Correlation or covariance matrix
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param modularity.hypot Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @param significance Logical. Indicates if goodness of fit test should be performed.
#' @param sample.size sample size in goodness of fit simulations via MonteCarlo 
#' @param iterations number os goodness of fit simulations
#' @param n number of jackknife resamplings
#' @param leave.out number of individuals to be left out of each jackknife, default is 10\%
#' @param ... aditional arguments to be passed to raply for the jackknife
#' @note Hypothesis can be named as column names, and these will be used to make labels in the 
#' output. 
#' 
#' Jackknife will return the best hypothesis for each sample.
#' @export
#' @aliases JackKnifeMINT
#' @references 
#' Marquez, E.J. 2008. A statistical framework for testing modularity in multidimensional data. 
#' Evolution 62:2688-2708. 
#' 
#' Parsons, K.J., Marquez, E.J., Albertson, R.C. 2012. Constraint and opportunity: the genetic 
#' basis and evolution of modularity in the cichlid mandible. The American Naturalist 179:64-78.
#' 
#' http://www-personal.umich.edu/~emarquez/morph/doc/mint_man.pdf
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
#' mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.1, 0.2) # between-modules
#' diag(mod.cor) = 1
#' mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
#' 
#' # True hypothesis and a bunch of random ones.
#' hypothetical.modules = cbind(modules, matrix(sample(c(1, 0), 4*15, replace=TRUE), 15, 4))
#' 
#' # if hypothesis columns are not named they are assigned numbers
#' colnames(hypothetical.modules) <- letters[1:7]
#'  
#' MINT(mod.cor, hypothetical.modules)
#' 
#' random_var = runif(15, 1, 10)
#' mod.data = mvtnorm::rmvnorm(100, sigma = sqrt(outer(random_var, random_var)) * mod.cor)
#' out_jack = JackKnifeMINT(mod.data, hypothetical.modules, n = 50)
#' 
#' library(ggplot2)
#' 
#' ggplot(out_jack, aes(rank, corrected.gamma)) + geom_point() + 
#'        geom_errorbar(aes(ymin = lower.corrected, ymax = upper.corrected))
MINT <- function(c.matrix, modularity.hypot, 
                            significance = FALSE, sample.size = NULL, iterations = 1000){
  cor.hypot.list = CombineHypot(modularity.hypot)
  if(is.list(cor.hypot.list)){
    if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- 1:dim(modularity.hypot)[2]
    gamma_df <- SingleMINT(cor.hypot.list, c.matrix, 
                               significance, sample.size, iterations)
    gamma_df = gamma_df[order(gamma_df$corrected.gamma),]
    rownames(gamma_df) <- NULL
    return(list(gamma_rank = gamma_df, 
                "Modularity_hypothesis" = modularity.hypot[,strsplit(gamma_df[1,1], "_")[[1]]]))
  } 
  else gamma(c.matrix, c.matrix * cor.hypot.list)
}

gamma <- function(S, S_0) sum(diag((x <- (S - S_0)) %*% t(x)))

SingleMINT <- function(cor.hypot.list, c.matrix, 
                           significance = FALSE, sample.size = NULL, iterations = 1000) {
  raw.gamma = laply(cor.hypot.list, function(x) gamma(c.matrix, c.matrix * x))
  num.zeros = laply(cor.hypot.list, function(x) sum(x[lower.tri(x)] == 0))  
  scaled.gamma = raw.gamma/max(raw.gamma)
  zeros_lm = lm(scaled.gamma ~ num.zeros)
  corrected.gamma = residuals(zeros_lm)
  gamma_df = data.frame(.id = names(cor.hypot.list), 
                        .n = 1:length(cor.hypot.list),
                        rank = rank(corrected.gamma),
                        corrected.gamma = corrected.gamma,
                        raw.gamma, stringsAsFactors = FALSE)
  if(significance){
    if(is.null(sample.size)) stop('sample.size must be non-null')
    gamma_df$sig = laply(cor.hypot.list, function(x){
      null_dist = MonteCarloStat(c.matrix * x, sample.size, 
                                 iterations, gamma, cov)$V1 < gamma(c.matrix * x, c.matrix)
      sum(null_dist)/iterations
    })
    # gamma_df$sig = laply(cor.hypot.list, function(x){
    #   tryCatch({null_dist = laply(rlply(iterations, solve(rwish(sample.size, solve(c.matrix * x))/sample.size)),
    #                               function(sim_x) gamma(sim_x, c.matrix * x)) < gamma(c.matrix * x, c.matrix)
    #   sum(null_dist)/iterations},
    #   error = function(c) { return(NA)})
    # })
  }
  return(gamma_df)
}

#' @export
#' @rdname MINT
JackKnifeMINT <- function(ind.data, modularity.hypot, 
                                     n = 1000, leave.out = floor(dim(ind.data)[1]/10),
                                     ...){
  if(isSymmetric(as.matrix(ind.data))) stop("input appears to be a matrix, use residuals")
  cor.hypot.list = CombineHypot(modularity.hypot)
  if(is.list(cor.hypot.list)){
    n.ind = dim(ind.data)[1]
    if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- 1:dim(modularity.hypot)[2]
    observed = SingleMINT(cor.hypot.list, cov(ind.data))
    out_jack = raply(n,function() {
      current.sample = sample(1:n.ind, n.ind - leave.out)
      as.matrix(SingleMINT(cor.hypot.list, cov(ind.data[current.sample,]))[,4:5])
    }, ...)
    jack_interval = aaply(out_jack, 2:3, quantile, c(0.025, 0.0975)) # jackknife confidence intervals
    out_order = apply(out_jack[,,1], 1, function(x) rank(x) == rank(observed$corrected.gamma))   # jackkkife ordering
    observed$rank.support = rowSums(out_order)/n
    observed$upper.raw = jack_interval[,2,2]
    observed$upper.corrected = jack_interval[,1,2]
    observed$lower.raw = jack_interval[,2,1]
    observed$lower.corrected = jack_interval[,1,1]
    observed = observed[order(observed$corrected.gamma),]
    rownames(observed) <- NULL
    return(observed)
  }
  else stop("use at least two hypothesis vectors, preferably more")
}
