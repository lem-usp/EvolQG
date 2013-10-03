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
