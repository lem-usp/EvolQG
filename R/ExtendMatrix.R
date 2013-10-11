#' Control Inverse matrix noise with Extension
#' 
#' Calculates the extented covariance matrix estimation
#' 
#' @param cov.matrix Covariance matrix
#' @param cut.off Cut off for second derivative variance. When below first below cut off, last eigen value is selected.
#' @return Extended covariance matrix.
#' @references Marroig, G., Melo, D. A. R., and Garcia, G. (2012). Modularity, noise, and natural selection. Evolution; international journal of organic evolution, 66(5), 1506-24. doi:10.1111/j.1558-5646.2011.01555.x
#' @author Diogo Melo
#' @note Covariance matrix being extended must be larger then 10x10
#' @examples 
#' cov.matrix = RandomMatrix(11, 1, 1, 100)
#' ext.matrix = ExtendMatrix(cov.matrix, cut.off = 10e-4)
#' @keyword extension
#' @keyword covariancematrix
ExtendMatrix <- function(cov.matrix, cut.off = NULL) {
  p = dim(cov.matrix)[1]
  if(p<10)
    stop("matrix is too small")
  eigen.cov.matrix = eigen(cov.matrix)
  eVal = eigen.cov.matrix$values
  eVec = eigen.cov.matrix$vectors
  if(is.null(cut.off)){
    grad = array(dim=c(p-2))
    tr.cov.matrix = sum(eVal)
    for (i in 1:(p-2))
      grad[i] = abs(eVal[i]/tr.cov.matrix - 2*(eVal[i+1]/tr.cov.matrix) + eVal[i+2]/tr.cov.matrix)
    var.grad = array(dim=c(p-6))
    for (i in 1:(p-6)){
      var.grad[i] = var(grad[i:(i+4)])
    }
    length(var.grad[var.grad<1e-4])
    x11()
    plot(4:(p-3),var.grad)
    cut.off = floor(locator(1)$x)
  }
  eVal[eVal < eVal[cut.off]] = eVal[cut.off]
  extended.cov.matrix = eVec%*%diag(eVal)%*%t(eVec)
  colnames(extended.cov.matrix) = colnames(cov.matrix)
  rownames(extended.cov.matrix) = rownames(cov.matrix)
  return(extended.cov.matrix)
}
