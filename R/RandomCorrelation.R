#' Random correlation matrix
#'
#' Internal function for generating random correlation matrices.
#' Use RandomMatrix() instead.
#' @param dimension Number of traits in random matrix
#' @param k Parameter for correlation matrix generation. Involves check for positive defitness
#' @return Random correlation matrix with size (dimension x dimension)
#' @author Edgar Zanella
#' @keywords randommatrices
#' @export
#' @useDynLib evolqg
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
RandomCorrelation <- function (dimension, k=10**-3) {
  m <- createRandomMatrix(dimension, k)
  return(m)
}
