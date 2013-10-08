RandCorr <- function(num.traits, ke = 10^-3){
  random.corr = array(0., dim = c(num.traits, num.traits))
  b = array(0., dim = c(num.traits, num.traits))
  b[lower.tri(b, diag = T)] = 1.
  random.corr[2:num.traits, 1] = -1 + 2*round(runif(num.traits-1)*10^8)/10^8
  b[2:num.traits, 1] = random.corr[2:num.traits, 1]
  for (i in 2:num.traits)
    b[i, 2:i] = sqrt(1 - random.corr[i, 1]^2)
  for (i in 3:num.traits){
    for (j in 2:(i-1)){
      b1 = b[j, 1:(j-1)]%*%b[i, 1:(j-1)]
      b2 = b[j, j]*b[i, j]
      z = b1 + b2
      y = b1 - b2
      if (b2 < ke){
        random.corr[i, j] = b1
        cosinv = 0
      }
      else
        random.corr[i, j] = y + (z - y)*round(runif(1)*10^8)/10^8
      cosinv = (random.corr[i, j] - b1)/b2
      if (is.finite(cosinv)){
        if (cosinv > 1)
          b[i, (j+1):num.traits] = 0
        else if (cosinv < -1){
          b[i, j] = -b[i, j]
          b[i, (j+1):num.traits] = 0
        }
        else{
          b[i, j] = b[i, j]*cosinv
          sinTheta = sqrt(1 - cosinv^2)
          for (k in (j+1):num.traits)
            b[i, k] = b[i, k]*sinTheta
        }
      }
    }
  }
  random.corr = random.corr + t(random.corr) + diag(rep(1, num.traits))
  perm = sample(1:num.traits)
  random.corr = (random.corr[perm,])[,perm]
  return (random.corr)
}

RandomMatrix <- function(num.traits, num.matrices = 1, min.var = 1, max.var = 1, variance = NULL, ke = 10^-3){
    if(num.matrices==1){
        if(is.null(variance)) variance <- runif(num.traits, min.var, max.var)
        rand.mat <- RandCorr(num.traits, ke) * sqrt(outer(variance, variance))
    }
    else{
        if(is.null(variance)) variance <- matrix(runif(num.matrices*num.traits, min.var, max.var), num.matrices, num.traits)
        else variance <-  matrix(rep(variance, each = num.matrices), num.matrices, num.traits)
        rand.mat <- lapply(as.list(1:num.matrices), function(x) RandCorr(num.traits, ke) * sqrt(outer(variance[x,], variance[x,])))
    }
    return(rand.mat)
}
