#' Mean Squared Correlations
#'
#' Calculates the mean squared correlation of a covariance or correlation matrix. Measures integration.
#'
#' @param c.matrix Covariance or correlation matrix.
#' @return Mean squared value of of diagonal elements of correlations
#' @export
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{Flexibility}}
#' @examples
#' cov.matrix <- RandomMatrix(10, 1, 1, 10)
#' CalcR2(cov.matrix)
#' CalcR2(cov2cor(cov.matrix))
#' @keywords correlation
#' @keywords integration

CalcR2 <- function (c.matrix){
    cor.matrix = cov2cor(c.matrix)
    return (mean (cor.matrix [lower.tri (cor.matrix)]^2))
}
