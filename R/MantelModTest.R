#' Test single modularity hypothesis using Mantel correlation
#'
#' Calculates the correlation and Mantel significance test between a hypothetical binary modularity
#' matrix and a correlation matrix. Also gives mean correlation within- and between-modules. 
#' This function is usually only called by TestModularity.
#' 
#' CalcAVG can be used when a significance test is not required.
#'
#' @param cor.hypothesis Hypothetical correlation matrix, with 1s within-modules and 0s between modules.
#' @param cor.matrix Observed empirical correlation matrix.
#' @param permutations Number of permutations used in significance calculation.
#' @param MHI Indicates if Modularity Hypothesis Indexshould be calculated instead of AVG Ratio.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param ... aditional arguments passed to MantelCor
#' @param landmark.dim Used if permutations should be performed mantaining landmark structure in geometric morphomotric data. Either 2 for 2d data or 3 for 3d data. Default is NULL for non geometric morphomotric data.
#' @param withinLandmark Logical. If TRUE within-landmark correlation are used in calculation of correlation. Only used if landmark.dim is passed, default is FALSE.
#' @return Returns a vector with the matrix correlation, significance via Mantel, within- and between module correlation.
#' @export
#' @importFrom vegan mantel
#' @importFrom Matrix bdiag
#' @rdname MantelModTest
#' @references Porto, Arthur, Felipe B. Oliveira, Leila T. Shirai, Valderes Conto, and Gabriel Marroig. 2009. "The Evolution of Modularity in the Mammalian Skull I: Morphological Integration Patterns and Magnitudes." Evolutionary Biology 36 (1): 118-35. doi:10.1007/s11692-008-9038-3.
#' 
#' Modularity and Morphometrics: Error Rates in Hypothesis Testing Guilherme Garcia, Felipe Bandoni de Oliveira, Gabriel Marroig bioRxiv 030874; doi: http://dx.doi.org/10.1101/030874
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{mantel}},\code{\link{MantelCor}},\code{\link{CalcAVG}},\code{\link{TestModularity}}
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
#' hypot.mask = matrix(as.logical(cor.hypot), 12, 12)
#' mod.cor = matrix(NA, 12, 12)
#' mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
#' mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
#' diag(mod.cor) = 1
#' mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
#' 
#' MantelModTest(cor.hypot, mod.cor)
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords manteltest
MantelModTest <- function (cor.hypothesis, cor.matrix, ...) UseMethod("MantelModTest")

#' @rdname MantelModTest
#' @method MantelModTest default
#' @export
MantelModTest.default <- function (cor.hypothesis, cor.matrix,
                                   permutations = 1000, MHI = FALSE, ..., 
                                   landmark.dim = NULL, withinLandmark = FALSE) {
  if(!any(cor.hypothesis == 0 | cor.hypothesis == 1)) stop("modularity hypothesis matrix should be binary")
  mantel.output <- MantelCor(cor.matrix, cor.hypothesis, permutations = permutations,
                             MHI, landmark.dim = landmark.dim, withinLandmark = FALSE, mod = TRUE, ...)
  output = c(mantel.output, 
             CalcAVG(cor.hypothesis, cor.matrix, MHI, landmark.dim))
  return (output)
}

#' @rdname MantelModTest
#' @method MantelModTest list
#' @export
MantelModTest.list <- function (cor.hypothesis, cor.matrix,
                                permutations = 1000, MHI = FALSE, 
                                landmark.dim = NULL, withinLandmark = FALSE,
                                ..., 
                                parallel = FALSE){
  output <- SingleComparisonMap(cor.hypothesis, cor.matrix,
                                function(x, y) MantelModTest(x, y,
                                                             permutations, MHI, 
                                                             landmark.dim = landmark.dim,
                                                             withinLandmark = withinLandmark),
                                parallel = parallel)
  return(output)
}
