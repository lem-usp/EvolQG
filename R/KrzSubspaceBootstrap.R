#' Quasi-Bayesian Krzanowski subspace comparison
#'
#' Calculates the usual Krzanowski subspace comparison using a posterior samples
#' for a set of phenotypic covariance matrices. Then, this observed comparison 
#' is contrasted to the subspace comparison across a permutation of the original
#' data. Residuals, which are used to calculate the observed P-matrices, are 
#' shuffled across groups. This process is repeated, creating a null distribution
#' of  subspace comparisons under the hypothesis that all P-matrices come from the
#' same population. This method is a modification on the fully Bayesian method 
#' proposed in Aguirre et. al 2014 and improved in Morrisey et al 2020.
#' 
#' @param x list of linear models from which P-matrices should be calculated
#' @param rep number of bootstrap samples to be made
#' @param MCMCsamples number of MCMCsamples for each P-matrix posterior distribution.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#'
#'
#'
#' @export
#' @examples 
#' 
#' 
#' 
#' 
KrzSubspaceBootstrap = function(x, rep = 1, MCMCsamples = 1000, parallel = FALSE){
  P_list = laply(x, function(x) BayesianCalculateMatrix(x, samples = MCMCsamples)$Ps)
  P_list = aperm(P_list, c(3, 4, 1, 2))
  res_list = lapply(x, residuals)
  n_list = sapply(res_list, nrow)

  residuals = do.call(rbind, res_list)

  Hs = llply(alply(P_list, 4, function(x) alply(x, 3)), function(x) KrzSubspace(x, 3)$H)
  avgH = Reduce("+", Hs)/length(Hs)
  avgH.vec <- eigen(avgH)$vectors
  MCMC.H.val = laply(Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))

  randomKrz = function(samples){
    random_lm = lapply(n_list, 
                       function(n) 
                         lm(residuals[sample(1:nrow(residuals), n),]~1))
    random_P_list = laply(random_lm, 
                          function(x) 
                             BayesianCalculateMatrix(x, samples = samples)$Ps)
    random_P_list = aperm(random_P_list, c(3, 4, 1, 2))
    Hs = llply(
      alply(random_P_list, 4, function(x) alply(x, 3)), 
      function(x) 
        KrzSubspace(x, 3)$H)
        avgH = Reduce("+", Hs)/length(Hs)
        avgH.vec <- eigen(avgH)$vectors
        MCMC.H.val.random = laply(Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))
        MCMC.H.val.random
  }

  rand = laply(1:rep, function(i) randomKrz(MCMCsamples), .parallel = parallel)

  MCMC.H.val.random = do.call(rbind, alply(rand, 1, identity))
  list(observed = MCMC.H.val, random = MCMC.H.val.random)
}

#' @importFrom coda HPDinterval
KrzSubspaceDataFrame <- function(observed, random, n = ncol(observed)){
  lower = NULL
  upper = NULL
  HPD.H.val <- cbind(HPDinterval(as.mcmc(observed)),
                     HPDinterval(as.mcmc(random)))

  dat.observed = as.data.frame(HPD.H.val[1:n,1:2])
  dat.random = as.data.frame(HPD.H.val[1:n,3:4])
  dat.observed$class = 'observed'
  dat.random$class = 'random'
  dat.observed$PC = dat.random$PC = paste('PC', 1:n, sep = '')
  dat = rbind(dat.observed, dat.random)
  dat$mean = rowMeans(cbind(lower, upper))

  orderlist = paste('PC', 1:n, sep = '')
  dat$PC = factor(dat$PC, orderlist)

  return(dat)
}

#' @importFrom ggplot2 ggplot geom_point geom_pointrange position_dodge aes xlab ylab
#' @export
PlotKrzSubspace = function(x){
  lower = NULL
  upper = NULL
  PC = NULL
  plot = ggplot(x, aes(PC, mean, color = class, group = interaction(PC, class))) +
    geom_point(position = position_dodge(width = 0.3)) +
    geom_pointrange(position = position_dodge(width = 0.3), aes(ymin=lower, ymax=upper)) +
    xlab('eigenvectors') + ylab('eigenvalues')
  return(plot)
}
