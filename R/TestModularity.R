#' Test modularity hypothesis
#'
#' Tests modularity hypothesis using cor.matrix matrix and trait groupings
#' @param cor.matrix Correlation matrix
#' @param modularity.hypot Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @param permutations Number of permutations, to be passed to MantelModTest
#' @param MHI Indicates if test should use Modularity Hypothesis Index instead of AVG Ratio
#' @param ... aditional arguments passed to MantelModTest
#' @param landmark.dim Used if permutations should be performed mantaining landmark structure in geometric morphomotric data. Either 2 for 2d data or 3 for 3d data. Default is NULL for non geometric morphomotric data.
#' @param withinLandmark Logical. If TRUE within-landmark correlations are used in the calculation of matrix correlation. Only used if landmark.dim is passed, default is FALSE.
#' @return Returns mantel correlation and associated probability for each modularity hypothesis, along with AVG+, AVG-, AVG Ratio for each module.
#' A total hypothesis combining all hypotesis is also tested.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MantelModTest}}
#' @export
#' @rdname TestModularity
#' @references Porto, Arthur, Felipe B. Oliveira, Leila T. Shirai, Valderes Conto, and Gabriel Marroig. 2009. "The Evolution of Modularity in the Mammalian Skull I: Morphological Integration Patterns and Magnitudes." Evolutionary Biology 36 (1): 118-35. doi:10.1007/s11692-008-9038-3.
#' @examples
#' cor.matrix <- RandomMatrix(10)
#' rand.hypots <- matrix(sample(c(1, 0), 30, replace=TRUE), 10, 3)
#' mod.test <- TestModularity(cor.matrix, rand.hypots)
#' 
#' cov.matrix <- RandomMatrix(10, 1, 1, 10)
#' cov.mod.test <- TestModularity(cov.matrix, rand.hypots, MHI = TRUE)
#' nosize.cov.mod.test <- TestModularity(RemoveSize(cov.matrix), rand.hypots, MHI = TRUE)
#' @keywords mantel
#' @keywords modularity
TestModularity <- function (cor.matrix, modularity.hypot, 
                            permutations = 1000, MHI = FALSE, ...,
                            landmark.dim = NULL, withinLandmark = FALSE) {
  if(is.null(dim(modularity.hypot))){
    return(MantelModTest(CreateHypotMatrix(modularity.hypot), cor.matrix, 
           permutations = permutations, MHI = MHI, 
           landmark.dim = landmark.dim, withinLandmark = withinLandmark, ...))
  }
  else{
    m.hyp.list <- CreateHypotMatrix(as.matrix(modularity.hypot))
    if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- 1:dim (modularity.hypot) [2]
    names(m.hyp.list) <- c(colnames(modularity.hypot),"Full Integration")
  }
  output <- MantelModTest(m.hyp.list, cor.matrix, 
                          permutations = permutations, MHI = MHI, 
                          landmark.dim = landmark.dim, withinLandmark = withinLandmark, ...)
  names(output)[1] <- 'hypothesis'
  return (output)
}