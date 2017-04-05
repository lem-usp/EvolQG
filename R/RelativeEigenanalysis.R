#' Relative Eigenanalysis
#'
#' Computes relative eigenvalues and eigenvectors between a pair of covariance matrices.
#'
#' @param cov.x covariance matrix
#' @param cov.y covariance matrix
#' @param symmetric compute symmetric eigenanalysis?
#' @return list with two objects: eigenvalues and eigenvectors
#'
#' @export
#' @importFrom expm sqrtm
#'
#' 
#' @rdname RelativeEigenanalysis
#'
#' @references Bookstein, F. L., and P. Mitteroecker, P. "Comparing Covariance Matrices by
#' Relative Eigenanalysis, with Applications to Organismal Biology." Evolutionary Biology 41, no. 2
#' (June 1, 2014): 336-350. doi:10.1007/s11692-013-9260-5.
#'
#' @author Guilherme Garcia, Diogo Melo
#'
#' @examples
#' data(dentus)
#' dentus.vcv <- dlply(dentus, .(species), function(df) var(df[, -5]))
#'
#' dentus.eigrel <- RelativeEigenanalysis(dentus.vcv [[1]], dentus.vcv[[5]])
#' 
RelativeEigenanalysis <- function(cov.x, cov.y, symmetric = FALSE)
{
    if(symmetric)
    {
        out.elem <- solve(sqrtm(cov.y))
        matrix.prod <- out.elem %*% cov.x %*% out.elem
    }
    else
        matrix.prod <- cov.x %*% solve(cov.y)

    eigen(matrix.prod)
}   
