#' Compare matrices via Mantel Correlation
#'
#' Calculates correlation matrix correlation and significance via Mantel test.
#'
#' @param cor.x Single correlation matrix or list of correlation matrices.
#'
#' If single matrix is suplied, it is compared to cor.y.
#'
#' If list is suplied and no cor.y is suplied, all matrices
#' are compared.
#'
#' If cor.y is suplied, all matrices in list are compared to it.
#' @param cor.y First argument is compared to cor.y.
#' Optional if cor.x is a list.
#' @param permutations Number of permutations used in significance calculation.
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @param ... aditional arguments passed to other methods
#' @param landmark.dim Used if permutations should be performed mantaining landmark structure in geometric morphomotric data. Either 2 for 2d data or 3 for 3d data. Default is NULL for non geometric morphomotric data.
#' @param withinLandmark Logical. If TRUE within-landmark correlations are used in the calculation of matrix correlation. Only used if landmark.dim is passed, default is FALSE.
#' @param mod Set TRUE to use mantel in testing modularity hypothesis. Should only be used in MantelModTest.
#' @return If cor.x and cor.y are passed, returns matrix pearson
#' correlation and significance via Mantel permutations.
#'
#' If cor.x is a list of matrices and cor.y is passed, same as above, but for all matrices in cor.x.
#'
#' If only cor.x is passed, a matrix of MantelCor average
#' values and probabilities of all comparisons.
#' If repeat.vector is passed, comparison matrix is corrected above
#' diagonal and repeatabilities returned in diagonal.
#' @note If the significance is not needed, MatrixCor provides the 
#' correlation and skips the permutations, so it is much faster.
#' @export
#' @importFrom vegan mantel
#' @rdname MantelCor
#' @references http://en.wikipedia.org/wiki/Mantel_test
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{KrzCor}},\code{\link{RandomSkewers}},\code{\link{mantel}},\code{\link{RandomSkewers}},\code{\link{TestModularity}}, \code{\link{MantelModTest}}
#' @examples
#' c1 <- RandomMatrix(10, 1, 1, 10)
#' c2 <- RandomMatrix(10, 1, 1, 10)
#' c3 <- RandomMatrix(10, 1, 1, 10)
#' MantelCor(cov2cor(c1), cov2cor(c2))
#' 
#' cov.list <- list(c1, c2, c3)
#' cor.list <- llply(list(c1, c2, c3), cov2cor)
#'
#' MantelCor(cor.list)
#'
#'# For repeatabilities we can use MatrixCor, which skips the significance calculation
#' reps <- unlist(lapply(cov.list, MonteCarloRep, 10, MatrixCor, correlation = TRUE))
#' MantelCor(cor.list, repeat.vector = reps)
#'
#' c4 <- RandomMatrix(10)
#' MantelCor(cor.list, c4)
#' 
#' #Multiple threads can be used with some foreach backend library, like doMC or doParallel
#' #library(doParallel)
#' ##Windows:
#' #cl <- makeCluster(2)
#' #registerDoParallel(cl)
#' ##Mac and Linux:
#' #registerDoParallel(cores = 2)
#' #MantelCor(cor.list, parallel = TRUE) 
#' @keywords matrixcomparison
#' @keywords matrixcorrelation
#' @keywords randomskewers
MantelCor <- function (cor.x, cor.y, ...) UseMethod("MantelCor")

#' @rdname MantelCor
#' @method MantelCor default
#' @export
MantelCor.default <- function (cor.x, cor.y, permutations = 1000, ..., 
                               landmark.dim = NULL, withinLandmark = FALSE, mod = FALSE) {
  if(!mod & (sum(diag(cor.x)) != dim(cor.x)[1] | sum(diag(cor.y))!= dim(cor.y)[1]))
    warning("Matrices do not appear to be correlation matrices. Use with caution.")
  if(!is.null(landmark.dim)){
    if(!any(landmark.dim == c(2, 3))) stop("landmark.dim should be either 2 or 3 dimensions")
    if(withinLandmark) cor.mask = lower.tri(cor.x, landmark.dim)
    else cor.mask = lower.tri.land(cor.x, landmark.dim)
    correlation <- cor(cor.x [cor.mask],
                       cor.y [cor.mask])
    null_vector = vector("numeric", permutations)
    num.traits = dim(cor.x)[1]
    n.land = num.traits/landmark.dim
    for(perm in 1:permutations){
      current.perm = rep(sample(1:n.land) * landmark.dim, each = landmark.dim) - (landmark.dim-1):0
      null_vector[perm] = cor(cor.x[cor.mask],
                              cor.y[current.perm, current.perm][cor.mask])
    }
    prob = sum(null_vector > correlation)/permutations
  } else{
    mantel.output <- mantel(cor.x, cor.y, permutations = permutations)
    correlation <- mantel.output$statistic
    prob <- mantel.output$signif
  }
  output <- c(correlation, prob)
  names(output) <- c("Rsquared", "Probability")
  return (output)
}

CreateWithinLandMat <- function(num.land, land.dim){
  num.traits = num.land * land.dim
  matrix(as.logical(bdiag(rlply(num.land, matrix(1, land.dim, land.dim)))), num.traits, num.traits)
}

lower.tri.land <- function(x, landmark.dim = NULL){
  num.land = dim(x)[1] / landmark.dim
  !CreateWithinLandMat(num.land, landmark.dim) & lower.tri(x)
}

#' @rdname MantelCor
#' @method MantelCor list
#' @export
MantelCor.list <- function (cor.x, cor.y = NULL,
                            permutations = 1000, repeat.vector = NULL,
                            parallel = FALSE, ...)
{
  if (is.null (cor.y)) {
    output <- ComparisonMap(cor.x,
                         function(x, cor.y) MantelCor(x, cor.y, permutations, ...),
                         repeat.vector = repeat.vector,
                         parallel = parallel)
  } else{
    output <- SingleComparisonMap(cor.x, cor.y,
                               function(x, y) MantelCor(y,
                                                        x,
                                                        permutations, ...),
                               parallel = parallel)
  }
  return(output)
}

#' @rdname MantelCor
#' @method MantelCor mcmc_sample
#' @export
MantelCor.mcmc_sample <- function (cor.x, cor.y, ..., parallel = FALSE)
{
  MatrixCor(cor.x, cor.y, parallel)
}

#' @export
#' @rdname MantelCor
MatrixCor <- function (cor.x, cor.y, ...) UseMethod("MatrixCor")

#' @rdname MantelCor
#' @method MatrixCor default
#' @export
MatrixCor.default <- function (cor.x, cor.y, ...)                           
{
  if(sum(diag(cor.x)) != dim(cor.x)[1] | sum(diag(cor.y))!= dim(cor.y)[1])
    warning("Matrices do not appear to be correlation matrices. Use with caution.")
  c("correlations" = cor(cor.x[lower.tri(cor.x)], cor.y[lower.tri(cor.y)], ...))
}

#' @rdname MantelCor
#' @method MatrixCor list
#' @export
MatrixCor.list <- function (cor.x, cor.y = NULL,
                            permutations = 1000, repeat.vector = NULL,
                            parallel = FALSE, ...)
{
  if (is.null (cor.y)) {
    output <- ComparisonMap(cor.x,
                            function(x, cor.y) c(MatrixCor(x, cor.y), NA),
                            repeat.vector = repeat.vector,
                            parallel = parallel)[[1]]
  } else{
    output <- SingleComparisonMap(cor.x, cor.y,
                                  function(x, y) MatrixCor(x, y),                                     
                                  parallel = parallel)
  }
  return(output)
}

#' @rdname MantelCor
#' @method MatrixCor mcmc_sample
#' @export
MatrixCor.mcmc_sample <- function (cor.x, cor.y, ..., parallel = FALSE)
{
  if (class (cor.y) == "mcmc_sample") {
    n = dim(cor.x)[1]
    if(dim(cor.y)[1] != n) stop("samples must be of same size")
    cor.x <- alply(cor.x, 1, cov2cor)
    output <- aaply(1:n, 1, function(i) MatrixCor(cor.x, 
                                                  cov2cor(cor.y[i,,]))$correlation,
                    .parallel = parallel)
    output <- as.numeric(output)
  } else{
    output <- SingleComparisonMap(alply(cor.x, 1), cor.y,
                                  function(x, y) MatrixCor(x, y),
                                  parallel = parallel)
  }
  return(output)
}