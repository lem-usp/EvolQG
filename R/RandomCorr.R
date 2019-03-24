#' Random correlation matrix
#'
#' Internal function for generating random correlation matrices.
#' Use RandomMatrix() instead.
#' @param num.traits Number of traits in random matrix
#' @param ke Parameter for correlation matrix generation. Involves check for positive defitness
#' @return Random Matrix
#' @author Diogo Melo Edgar Zanella
#' @importFrom stats rbeta
#' @keywords randommatrices
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
      } else
        random.corr[i, j] = y + (z - y)*round(runif(1)*10^8)/10^8
      cosinv = c((random.corr[i, j] - b1)/b2)
      if (is.finite(cosinv)){
        if (cosinv > 1)
          b[i, (j+1):num.traits] = 0
        else if (cosinv < -1){
          b[i, j] = -b[i, j]
          b[i, (j+1):num.traits] = 0
        } else {
          b[i, j] = b[i, j]*cosinv
          sinTheta = sqrt(1 - cosinv^2)

          b[i, (j+1):num.traits] = b[i, (j+1):num.traits]* c(sinTheta)
        }
      }
    }
  }
  random.corr = random.corr + t(random.corr) + diag(rep(1, num.traits))
  perm = sample(1:num.traits)
  random.corr = (random.corr[perm,])[,perm]
  return (random.corr)
}

# https://stats.stackexchange.com/questions/2746/how-to-efficiently-generate-random-positive-semidefinite-correlation-matrices/125017
RandLKJ = function(num.traits, shape){
  P = matrix(0, num.traits, num.traits)       
  random.corr = diag(num.traits)

  for (k in 1:(num.traits-1)){
    for (i in (k+1):num.traits){
      P[k,i] = rbeta(1, shape, shape) # sampling from beta
      P[k,i] = (P[k,i]-0.5)*2         # linearly shifting to [-1, 1]
      p = P[k,i]
      for (l in seq((k-1), 1)){       # converting partial correlation to raw correlation
        if(l != 0)
          p = p * sqrt((1-P[l,i]^2) * (1-P[l,k]^2)) + P[l,i]*P[l,k]
      }
      random.corr[k,i] = p
      random.corr[i,k] = p
    }
  }
  
  perm = sample(1:num.traits)
  random.corr = (random.corr[perm,])[,perm]
  
  return(random.corr)
}
