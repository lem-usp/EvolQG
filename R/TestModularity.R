#' Test modularity hypothesis
#'
#' Tests modularity hypothesis using cor.matrix matrix and trait groupings
#' @param cor.matrix Correlation matrix
#' @param modularity.hypot Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @param permutations Number of permutations, to be passed to MantelModTest
#' @param MHI Indicates if test should use Modularity Hypothesis Index instead of AVG Ratio
#' @param ... aditional arguments passed to MantelModTest
#' @param landmark.dim Used if permutations should be performed mantaining landmark structure in geometric morphomotric data. Either 2 for 2d data or 3 for 3d data. Default is NULL for non-geometric morphomotric data.
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
TestModularity <- function (cor.matrix, modularity.hypot, permutations = 1000, MHI = FALSE, ...,
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
  output <- MantelModTest(m.hyp.list, cor.matrix, permutations = permutations, MHI = MHI, 
                          landmark.dim = landmark.dim, withinLandmark = withinLandmark, ...)
  names(output)[1] <- 'hypothesis'
  return (output)
}

#' @export
#' @rdname TestModularity
CreateHypotMatrix <- function(modularity.hypot){
  if(is.null(dim(modularity.hypot))) return(outer(modularity.hypot, modularity.hypot))
  num.hyp <- dim (modularity.hypot) [2]
  num.traits <- dim (modularity.hypot) [1]
  m.hyp.list <- alply(modularity.hypot, 2, function(x) outer(x, x))
  m.hyp.list[[num.hyp+1]] <- matrix(as.integer (as.logical (Reduce ("+", m.hyp.list[1:num.hyp]))),
                                    num.traits, num.traits, byrow=T)
  return(m.hyp.list[1:(num.hyp+1)])
}

#' Test single modularity hypothesis using Mantel correlation
#'
#' Calculates the correlation and Mantel significance test between a hypothetical binary modularity
#' matrix and a correlation matrix. Also gives mean correlation within- and between-modules. 
#' This function is usually only called by TestModularity.
#' 
#' CalcAVG can be used when a significance test is not required.
#'
#' @param cor.hypothesis hypothetical correlation matrix, with 1s within-modules and 0s between modules
#' @param cor.matrix Observed empirical correlation matrix.
#' @param permutations Number of permutations used in significance calculation.
#' @param MHI Indicates if modularity hypothesis test should use Modularity Hypothesis Index instead of AVG Ratio.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param ... aditional arguments passed to MantelCor
#' @param landmark.dim Used permutations should be performed mantaining landmark structure in geometric morphomotric data. Either 2 for 2d data or 3 for 3d data. Default is NULL for non-geometric morphomotric data.
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
#' @seealso \code{\link{mantel}},\code{\link{MantelCor}},\code{\link{TestModularity}}
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

#' @rdname MantelModTest
#' @export
CalcAVG <- function(cor.hypothesis, cor.matrix, MHI = FALSE, landmark.dim = NULL){
  if(!all(cor.hypothesis == 0 | cor.hypothesis == 1)) stop("modularity hypothesis matrix should be binary")
  if(!is.null(landmark.dim)){
    if(!any(landmark.dim == c(2, 3))) stop("landmark.dim should be either 2 or 3 dimensions")
    num.traits = dim(cor.matrix)[1]
    n.land = num.traits/landmark.dim
    withinLandMat = CreateWithinLandMat(n.land, landmark.dim)
    cor.hypothesis[withinLandMat] = 2
  }
  index <- cor.hypothesis[lower.tri(cor.hypothesis)]
  avg.plus <- mean (cor.matrix [lower.tri(cor.matrix)] [index == 1])
  avg.minus <- mean (cor.matrix [lower.tri(cor.matrix)] [index == 0])
  if(MHI){
    avg.index <- (avg.plus - avg.minus)/CalcICV(cor.matrix)
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


gamma <- function(S, S_0) sum(diag((x <- (S - S_0)) %*% t(x)))

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
#' ModularityGamma(mod.cor, hypot)
ModularityGamma <- function(cor.matrix, hypot){
  cor.hypot = CreateHypotMatrix(hypot)
  if(class(cor.hypot) == "list") cor.hypot = cor.hypot[[length(cor.hypot)]]
  diag(cor.hypot) <- 1
  gamma(cor.matrix, cor.matrix * cor.hypot)
}

#'@examples
#' rand.hypots <- matrix(sample(c(1, 0), 20, replace=TRUE), 10, 2)
#' CombineHypot(rand.hypots)
CombineHypot <- function(modularity.hypot){
  n.hypots = dim(modularity.hypot)[2]  
  counter = BinToDec(rep(1, n.hypots))
  hypot_list = list()
  k = 1
  for(i in seq(counter)){
    mask = DecToBin(i)
    mask = as.logical(as.numeric((mask[(32-(n.hypots-1)):32])))
    if(sum(mask) > 1) new_hypot = CreateHypotMatrix(modularity.hypot[,mask])[[sum(mask)+1]]
    else new_hypot = CreateHypotMatrix(modularity.hypot[,mask])
    diag(new_hypot) <- 1
    if(length(hypot_list) == 0){
      hypot_list[[k]] = new_hypot
      k = k + 1
    }
    else if(!any(laply(hypot_list, function(x) all(x == new_hypot)))){ 
      hypot_list[[k]] = new_hypot
      k = k + 1
    }
  }
  hypot_list
}

# http://stackoverflow.com/questions/12892348/convert-binary-string-to-binary-or-decimal-value
BinToDec <- function(x) sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
# http://stackoverflow.com/questions/6614283/converting-decimal-to-binary-in-r
DecToBin <- function(x) sapply(strsplit(paste(rev(intToBits(x))),""),`[[`,2)

