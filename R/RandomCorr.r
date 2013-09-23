RandCorr  <- function(n, ke = 10^-3){
  random.corr = array(0., dim = c(n, n))
  b = array(0., dim = c(n, n))
  b[lower.tri(b, diag = T)] = 1.
  random.corr[2:n, 1] = -1 + 2*round(runif(n-1)*10^8)/10^8
  b[2:n, 1] = random.corr[2:n, 1]
  for (i in 2:n)
    b[i, 2:i] = sqrt(1 - random.corr[i, 1]^2)
  for (i in 3:n){
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
          b[i, (j+1):n] = 0
        else if (cosinv < -1){
          b[i, j] = -b[i, j]
          b[i, (j+1):n] = 0
        }
        else{
          b[i, j] = b[i, j]*cosinv
          sinTheta = sqrt(1 - cosinv^2)
          for (k in (j+1):n)
            b[i, k] = b[i, k]*sinTheta
        }
      }
    }
  }
  random.corr = random.corr + t(random.corr) + diag(rep(1, n))
  perm = sample(1:n)
  random.corr = (random.corr[perm,])[,perm]
  return (random.corr)
}
