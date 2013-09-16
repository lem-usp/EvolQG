rmvNorm2 <-
function (n, theta = rep(0, nrow(sigma)),
                      sigma = diag(length(theta)),
                      method = c("eigen", "svd", "chol"))
  # Calculates random deviates from a Normal multivariate distribution
  #
  # Args:
  #   n: number os deviates
  #   theta: vetor of means
  #   sigma: covariance matrix
  #   method: generation method
  #
  # Return:
  #   Vector of deviates
{
  if (length(theta) != nrow(sigma)) {
    stop("theta and sigma have non-conforming size")
  }
  sigma.aux <- sigma
  dimnames(sigma.aux) <- NULL
  if (!isTRUE(all.equal(sigma.aux, t(sigma.aux)))) {
    warning("sigma is numerically not symmetric")
  }
  method <- match.arg(method) #  what is this?
  if (method == "eigen") {
    sigma.eigen <- eigen(sigma, symmetric = TRUE)
    if (!all(sigma.eigen$values >=
             -sqrt(.Machine$double.eps) * abs(sigma.eigen$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    random.deviates <- sigma.eigen$vectors %*%
                       diag(sqrt(sigma.eigen$values),
                            length(sigma.eigen$values)) %*%
                       t(sigma.eigen$vectors)
  }
  else if (method == "svd") {
    sigma.svd <- svd(sigma)
    if (!all(sigma.svd$d >= -sqrt(.Machine$double.eps) * abs(sigma.svd$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    random.deviates <- t(sigma.svd$v %*% (t(sigma.svd$u) * sqrt(sigma.svd$d)))
  }
  else if (method == "chol") {
    sigma.chol <- chol(sigma)
    random.deviates <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% sigma.chol
  }
  random.deviates <- sweep(random.deviates, 2, theta, "+")
  return(random.deviates)
}
