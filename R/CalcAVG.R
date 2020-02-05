#' Calculates mean correlations within- and between-modules
#'
#' Uses a binary correlation matrix as a mask to calculate average whitin- and between-module
#' correlations. Also calculates the ratio between them and the Modularity Hypothesis Index.
#' @param cor.hypothesis Hypothetical correlation matrix, with 1s within-modules and 0s between modules
#' @param cor.matrix Observed empirical correlation matrix.
#' @param MHI Indicates if Modularity Hypothesis Index should be calculated instead of AVG Ratio.
#' @param landmark.dim Used if within-landmark correlations are to be excluded in geometric morphomotric data. Either 2 for 2d data or 3 for 3d data. Default is NULL for non geometric morphomotric data.
#' @return a named vector with the mean correlations and derived statistics
#' @export
#' @examples
#' # Module vectors
#' modules = matrix(c(rep(c(1, 0, 0), each = 5),
#'                    rep(c(0, 1, 0), each = 5),
#'                    rep(c(0, 0, 1), each = 5)), 15)
#'
#' # Binary modular matrix
#' cor.hypot = CreateHypotMatrix(modules)[[4]]
#'
#' # Modular correlation matrix
#' hypot.mask = matrix(as.logical(cor.hypot), 15, 15)
#' mod.cor = matrix(NA, 15, 15)
#' mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
#' mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
#' diag(mod.cor) = 1
#' mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
#'
#' CalcAVG(cor.hypot, mod.cor)
#' CalcAVG(cor.hypot, mod.cor, MHI = TRUE)
CalcAVG <- function(cor.hypothesis, cor.matrix, MHI = TRUE, landmark.dim = NULL){
  if(!all(cor.hypothesis == 0 | cor.hypothesis == 1)) stop("modularity hypothesis matrix should be binary")
  if(!is.null(landmark.dim)){
    if(!any(landmark.dim == c(2, 3))) stop("landmark.dim should be either 2 or 3 dimensions")
    num.traits = dim(cor.matrix)[1]
    n.land = num.traits/landmark.dim
    withinLandMat = CreateWithinLandMat(n.land, landmark.dim)
    cor.hypothesis[withinLandMat] = 2
  }
  index <- cor.hypothesis[lower.tri(cor.hypothesis)]
  avg.plus  <- mean (cor.matrix [lower.tri(cor.matrix)] [index == 1])
  avg.minus <- mean (cor.matrix [lower.tri(cor.matrix)] [index == 0])
  if(MHI){
    avg.index <- (avg.plus - avg.minus)/CalcEigenVar(cor.matrix,sd = TRUE,rel=FALSE)
    output <- c(avg.plus, avg.minus, avg.index)
    names(output) <- c("AVG+", "AVG-", "MHI")
  } else{
    avg.ratio <- avg.plus / avg.minus
    output <- c(avg.plus, avg.minus, avg.ratio)
    names(output) <- c("AVG+", "AVG-", "AVG Ratio")
  }
  if(!is.null(landmark.dim)){
    output = c(output, mean (cor.matrix [lower.tri(cor.matrix)] [index == 2]))
    names(output)[length(output)] = "AVG within landmark"
  }
  return(output)
}