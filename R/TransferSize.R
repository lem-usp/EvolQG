TransferSize <-
function (cov.matrix.1, cov.matrix.2)
  # Transfers the size variation between matrices
  #
  # Args:
  #   cov.matrix{1,2}: A simetric covariance matrix
  # Return:
  #   A list with the matrices with the sizes swiched
{
  eigen.cov.matrix.1 = eigen(cov.matrix.1)
  eVal.1 = eigen.cov.matrix.1$values
  eVec.1 = eigen.cov.matrix.1$vectors

  eigen.cov.matrix.2 = eigen(cov.matrix.2)
  eVal.2 = eigen.cov.matrix.2$values
  eVec.2 = eigen.cov.matrix.2$vectors

  aux = eVal.1[1]
  eVal.1[1] = eVal.2[1]
  eVal.2[1] = aux

  cov.matrix.1.size.transfer <- eVec.1%*%diag(eVal.1)%*%t(eVec.1)
  cov.matrix.2.size.transfer <- eVec.2%*%diag(eVal.2)%*%t(eVec.2)

  return (list(cov.matrix.1.size.transfer, cov.matrix.2.size.transfer))
}
