KrzProjection <- function(cov.matrix, ...) UseMethod("KrzProjection")

KrzProjection.default <- function (cov.matrix.1, cov.matrix.2,
                           ret.dim.1 = NULL, ret.dim.2 = NULL)
  # Calculates the modified Krzanowski correlation between matrices,
  # projecting the variance in each principal  components of the first
  # matrix in to the ret.dim.2 components of the second.
  #
  # Args:
  #     cov.matrix.(1,2): covariance being compared
  #     ret.dim.2: number of retained dimensions in the comparison,
  #              default for nxn matrix is n/2-1
  # Return:
  #     Ratio of projected variance to total variance, and ratio of projectecd total in each PC
{
  require(plyr)
  num.traits <- dim(cov.matrix.1)[1]
  if (is.null(ret.dim.1))
    ret.dim.1 <- num.traits/2 - 1
  if (is.null(ret.dim.2))
    ret.dim.2 <- num.traits/2 - 1
  eigen.cov.matrix.1 <- eigen(cov.matrix.1)
  eVal.1 <- eigen.cov.matrix.1$values
  eVec.1 <- eigen.cov.matrix.1$vectors
  eVar.1 <- t(aaply(1:num.traits, 1, function(n) eVec.1[,n]*sqrt(eVal.1[n])))
  eVec.2 <- eigen(cov.matrix.2)$vectors
  SumSq  <- function(x) sum(x^2)
  MapProjection  <- function(x) SumSq(t(aaply(1:ret.dim.2, 1, function(n) eVar.1[,x]%*%eVec.2[,n])))
  ProjectionNorms <- aaply(1:ret.dim.1, 1, MapProjection)
  output <- list(total.variation = sum(ProjectionNorms)/sum(eVal.1),
                 per.PC = ProjectionNorms/eVal.1[1:ret.dim.1])
  return (output)
}

KrzProjection.list <- function(matrix.list,
                               ret.dim.1 = NULL, ret.dim.2 = NULL,
                               num.cores = 1,
                               full.results = FALSE){
  # Calculates the modified Krzanowski correlation between matrices,
  # projecting the variance in each principal  components of the first
  # matrix in to the ret.dim.2 components of the second for a list of matrices
  #
  # Args:
  #     cov.matrix.(1,2): covariance being compared
  #     ret.dim.2: number of retained dimensions in the comparison,
  #              default for nxn matrix is n/2-1
  #     num.cores: Number of cores for parallel computation. Requires doMC and foreach
  #     full.results: if FALSE returns only total variance, if TRUE also per PC variance.
  # Return:
  #     Ratio of projected variance to total variance, and ratio of projectecd total in each PC
  require(plyr)
  require(reshape2)
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else
    parallel = FALSE
  if(full.results){
    CompareToNProj <- function(n) llply(matrix.list,
                                        function(x) {KrzProjection(x,
                                                                   matrix.list[[n]],
                                                                   ret.dim.1, ret.dim.2)},
                                        .parallel = parallel)
  }
  else{
    CompareToNProj <- function(n) llply(matrix.list[names(matrix.list) != n],
                                        function(x) {KrzProjection(x,
                                                                   matrix.list[[n]],
                                                                   ret.dim.1, ret.dim.2)[1]},
                                        .parallel = parallel)
  }
  if(is.null(names(matrix.list))) {names(matrix.list) <- 1:length(matrix.list)}
  comparisons.proj <- llply(names(matrix.list),
                             CompareToNProj,
                             .parallel = parallel)
  if(full.results){
    names(comparisons.proj) = names(matrix.list)
    return(comparisons.proj)
  }
  else{
    comparisons.proj <- melt(comparisons.proj)
    comparisons.proj[,4] = names(matrix.list)[(comparisons.proj[,4])]
    comparisons.proj = comparisons.proj[,-2]
    comparisons.proj = acast(comparisons.proj, L2~L1)[names(matrix.list), names(matrix.list)]
    diag(comparisons.proj) = 0.
  }
  return(comparisons.proj)
}
