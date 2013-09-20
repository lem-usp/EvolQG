MultiKrzProjection <-
function(mat.list,
                               ret.dim.1 = NULL, ret.dim.2 = NULL,
                               num.cores = 1,
                               full.results = F){
  # Calculates the modified Krzanowski correlation between matrices,
  # projecting the variance in each principal  components of the first
  # matrix in to the ret.dim.2 components of the second.for a list of matrices
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
    CompareToNProj <- function(n) llply(mat.list,
                                        function(x) {KrzProjection(x,
                                                                   mat.list[[n]],
                                                                   ret.dim.1, ret.dim.2)},
                                        .parallel = parallel)
  }
  else{
    CompareToNProj <- function(n) llply(mat.list[names(mat.list) != n],
                                        function(x) {KrzProjection(x,
                                                                   mat.list[[n]],
                                                                   ret.dim.1, ret.dim.2)[1]},
                                        .parallel = parallel)
  }
  if(is.null(names(mat.list))) {names(mat.list) <- 1:length(mat.list)}
  comparisons.proj <- llply(names(mat.list),
                             CompareToNProj,
                             .parallel = parallel)
  if(full.results){
    names(comparisons.proj) = names(mat.list)
    return(comparisons.proj)
  }
  else{
    comparisons.proj <- melt(comparisons.proj)
    comparisons.proj[,4] = names(mat.list)[(comparisons.proj[,4])]
  }
  return(comparisons.proj[,-2])
}
