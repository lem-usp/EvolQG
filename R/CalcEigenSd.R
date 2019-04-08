#' Integration
#'
#' Calculates the standardized variance of the eigenvalues of a correlation matrix. 
#' This is a measure of overall integration between traits.
#'
#' @param c.matrix Covariance or correlation matrix.
#' @return Mean squared value of off diagonal elements of correlation matrix
#' @export
#' @author Diogo Melo
#' @seealso \code{\link{CalcR2}}, \code{\link{CalcICV}}, \code{\link{Flexibility}}
#' @references Pavlicev, Mihaela, James M. Cheverud, and Gunter P. Wagner. 2009. "Measuring Morphological Integration Using Eigenvalue Variance." Evolutionary Biology 36 (1): 157-70.
#' @examples
#' cov.matrix <- RandomMatrix(10, 1, 1, 10)
#' # both of the following calls are equivalent,
#' # CalcEigenSd() converts covariance matrices to correlation matrices internally
#' CalcEigenSd(cov.matrix)
#' CalcEigenSd(cov2cor(cov.matrix))
#' @keywords correlation
#' @keywords integration
CalcEigenSd <- function(c.matrix){
  cor.matrix = cov2cor(c.matrix)
  p = nrow(cor.matrix)
  eVals = eigen(cor.matrix)$values
  var_eVals = var(eVals) * (p - 1) / p
  return (sqrt(var_eVals)/sqrt(p-1))
}