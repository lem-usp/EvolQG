KrzProjection <- function(cov.x, cov.y, ...) UseMethod("KrzProjection")

KrzProjection.default <- function (cov.x, cov.y, ret.dim.1 = NULL, ret.dim.2 = NULL, ...) {
  require(plyr)
  num.traits <- dim(cov.x)[1]
  if (is.null(ret.dim.1))
    ret.dim.1 <- num.traits/2 - 1
  if (is.null(ret.dim.2))
    ret.dim.2 <- num.traits/2 - 1
  eigen.cov.x <- eigen(cov.x)
  eVal.1 <- eigen.cov.x$values
  eVec.1 <- eigen.cov.x$vectors
  eVar.1 <- t(aaply(1:num.traits, 1, function(n) eVec.1[,n]*sqrt(eVal.1[n])))
  eVec.2 <- eigen(cov.y)$vectors
  SumSq  <- function(x) sum(x^2)
  MapProjection  <- function(x) SumSq(t(aaply(1:ret.dim.2, 1, function(n) eVar.1[,x]%*%eVec.2[,n])))
  ProjectionNorms <- aaply(1:ret.dim.1, 1, MapProjection)
  output <- list(total.variation = sum(ProjectionNorms)/sum(eVal.1),
                 per.PC = ProjectionNorms/eVal.1[1:ret.dim.1])
  return (output)
}

KrzProjection.list <- function(cov.x, cov.y = NULL,
                               ret.dim.1 = NULL, ret.dim.2 = NULL,
                               num.cores = 1, full.results = FALSE, ...){
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
    CompareToNProj <- function(n) llply(cov.x,
                                        function(x) {KrzProjection(x,
                                                                   cov.x[[n]],
                                                                   ret.dim.1, ret.dim.2)},
                                        .parallel = parallel)
  }
  else{
    CompareToNProj <- function(n) llply(cov.x[names(cov.x) != n],
                                        function(x) {KrzProjection(x,
                                                                   cov.x[[n]],
                                                                   ret.dim.1, ret.dim.2)[1]},
                                        .parallel = parallel)
  }
  if(is.null(names(cov.x))) {names(cov.x) <- 1:length(cov.x)}
  comparisons.proj <- llply(names(cov.x),
                             CompareToNProj,
                             .parallel = parallel)
  if(full.results){
    names(comparisons.proj) = names(cov.x)
    return(comparisons.proj)
  }
  else{
    comparisons.proj <- melt(comparisons.proj)
    comparisons.proj[,4] = names(cov.x)[(comparisons.proj[,4])]
    comparisons.proj = comparisons.proj[,-2]
    comparisons.proj = acast(comparisons.proj, L2~L1)[names(cov.x), names(cov.x)]
    diag(comparisons.proj) = 0.
  }
  return(comparisons.proj)
}
