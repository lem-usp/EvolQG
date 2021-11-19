frobenius.norm <- function(x){
  return(entrywise.norm(x, 2))
}

entrywise.norm <- function(x, p) {
    if (!is.numeric(x)) {
      stop("argument x is not numeric")
    }
    if (is.matrix(x)) {
      Xmat <- x
    }
    else {
      if (is.vector(x)) {
        Xmat <- matrix(x, nrow = length(x), ncol = 1)
      }
      else {
        stop("argument x is neither vector nor matrix")
      }
    }
    if (p == 0) {
      stop("exponent p is zero")
    }
    if (is.infinite(p)) {
      return(maximum.norm(x))
    }
    return((sum(abs(Xmat)^p))^(1/p))
  }

maximum.norm <- function (x) {
  if (!is.numeric(x)) {
    stop("argument x is not numeric")
  }
  if (is.matrix(x)) {
    Xmat <- x
  }
  else {
    if (is.vector(x)) {
      X.mat <- x
    }
    else {
      stop("argument is neither a matrix nor a vector")
    }
  }
  norm <- max(abs(Xmat))
  return(norm)
}

frobenius.prod <- function (x, y){
  return(sum(hadamard.prod(x, y)))
}

hadamard.prod <- function (x, y) 
{
  if (!is.numeric(x)) {
    stop("argument x is not numeric")
  }
  if (!is.numeric(y)) {
    stop("argument y is not numeric")
  }
  if (is.matrix(x)) {
    Xmat <- x
  }
  else {
    if (is.vector(x)) {
      Xmat <- matrix(x, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (is.matrix(y)) {
    Ymat <- y
  }
  else {
    if (is.vector(y)) {
      Ymat <- matrix(y, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (nrow(Xmat) != nrow(Ymat)) 
    stop("argumentx x and y do not have the same row order")
  if (ncol(Xmat) != ncol(Ymat)) 
    stop("arguments x and y do not have the same column order")
  return(Xmat * Ymat)
}