Autonomy <- function (beta, cov.matrix) return ((1/(t (beta) %*% solve (cov.matrix, beta))) / (t (beta) %*% cov.matrix %*% beta))
ConditionalEvolvability <- function (beta, cov.matrix) return (1/(t (beta) %*% solve (cov.matrix, beta)))
Constraints <- function (beta, cov.matrix) return (abs (t (Normalize (eigen (cov.matrix)$vectors[,1])) %*% Normalize (cov.matrix %*% beta)))
Evolvability <- function (beta, cov.matrix) return (t (beta) %*% cov.matrix %*% beta)
Flexibility <- function (beta, cov.matrix) return (t (beta) %*% cov.matrix %*% beta / Norm (cov.matrix %*% beta))
Pc1Percent <- function (cov.matrix) return (eigen (cov.matrix)$values [1] / sum (eigen (cov.matrix)$values))
Respondability <- function (beta, cov.matrix) return (Norm (cov.matrix %*% beta))

MeanMatrixStatistics <- function (cov.matrix, iterations = 1000, full.results = F, num.cores = 1) {
  library(plyr)
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else{
    parallel = FALSE
  }
  matrix.stat.functions = list ('respondability' = Respondability,
                                'evolvability' = Evolvability,
                                'conditional.evolvability' = ConditionalEvolvability,
                                'autonomy' = Autonomy,
                                'flexibility' = Flexibility,
                                'constraints' = Constraints)
  num.traits <- dim (cov.matrix) [1]
  beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
  beta.mat <- apply (beta.mat, 2, Normalize)
  iso.vec <- Normalize (rep(1, num.traits))
  null.dist <- abs (t (iso.vec) %*% beta.mat)
  null.dist <- sort (null.dist)
  crit.value <- null.dist [round (0.95 * iterations)]
  cat ('critical value: ', crit.value, '\n')
  MatrixStatisticsMap <- function (CurrentFunc) return (apply (beta.mat, 2, CurrentFunc, cov.matrix = cov.matrix))
  stat.dist <- t(laply (matrix.stat.functions, MatrixStatisticsMap, .parallel = parallel))
  stat.dist <- cbind (stat.dist, null.dist)
  colnames (stat.dist) <- c('respondability',
                            'evolvability',
                            'conditional.evolvability',
                            'autonomy',
                            'flexibility',
                            'constraints',
                            'null.dist')
  stat.mean <- colMeans (stat.dist[,-7])
  integration <- c (CalcR2 (cov.matrix), Pc1Percent (cov.matrix))
  names (integration) <- c ('MeanSquaredCorrelation', 'pc1%')
  stat.mean <- c (integration, stat.mean)
  if(full.results)
    return (list ('dist' = stat.dist, 'mean' = stat.mean))
  else
    return (stat.mean)
}
