KrzProjection <-
function (cov.matrix.1, cov.matrix.2,
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
