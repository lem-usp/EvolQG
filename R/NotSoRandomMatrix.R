NotSoRandomMatrix <-
function (size, first.evec = NULL, eigen.values = NULL)
  {
    # generates a random covariance matrix
    #
    # args:
    #   size: number of traits
    #   first.evec: first eigenvector (optional)
    #   eigen.values: eigenvalue distribution (optional)
    # returns:
    #   random covariance matrix (or quite so, if extra arguments are provided)
    ei.vecs <- array (rnorm (size ^ 2), c(size, size))
    if (!is.null (first.evec))
      ei.vecs [,1] = first.evec
    ei.vecs <- apply (ei.vecs, 2, Normalize)
    if (!is.null (eigen.values))
      ei.vals <- sort (eigen.values, TRUE)
    else
      {
        ei.vals <- sort (rchisq (size, 1, 0), TRUE)
        ei.vals <- ei.vals / mean (ei.vals)
      }
    ei.vecs <- qr.Q (qr (ei.vecs))
    rand.mat <- ei.vecs %*% diag (ei.vals) %*% t (ei.vecs)
    return (rand.mat)
  }
