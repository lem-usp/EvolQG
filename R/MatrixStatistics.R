Autonomy <- function (beta, cov.matrix) return (((t (beta) %*% solve (cov.matrix, beta))^(-1)) / (t (beta) %*% cov.matrix %*% beta))
ConditionalEvolvability <- function (beta, cov.matrix) return ((t (beta) %*% solve (cov.matrix) %*% beta)^(-1))
Constraints <- function (beta, cov.matrix) return (abs (t (Normalize (eigen (cov.matrix)$vectors[,1])) %*% Normalize (cov.matrix %*% beta)))
Evolvability <- function (beta, cov.matrix) return (t (beta) %*% cov.matrix %*% beta)
Flexibility <- function (beta, cov.matrix) return (t (beta) %*% cov.matrix %*% beta / Norm (cov.matrix %*% beta))
Pc1Percent <- function (cov.matrix) return (eigen (cov.matrix)$values [1] / sum (eigen (cov.matrix)$values))
Respondability <- function (beta, cov.matrix) return (Norm (cov.matrix %*% beta))

MeanMatrixStatistics <- function (cov.matrix, iterations = 10000) {
  matrix.stat.functions = list ('respondability' = Respondability,
                                'evolvability' = Evolvability,
                                'conditional.evolvability' = ConditionalEvolvability,
                                'autonomy' = Autonomy,
                                'flexibility' = Flexibility,
                                'constraints' = Constraints)
  num.traits <- dim (cov.matrix) [1]
  beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
  beta.mat <- apply (beta.mat, 2, Normalize)
  iso.vec <- Normalize (rep(1, times = num.traits))
  null.dist <- abs (t (iso.vec) %*% beta.mat)
  null.dist <- sort (null.dist)
  crit.value <- null.dist [round (0.95 * iterations)]
  cat ('critical value: ', crit.value, '\n')
  stat.dist <- array (0, c(iterations, 8))
  MatrixStatisticsMap <- function (CurrentFunc) return (apply (beta.mat, 2, CurrentFunc, cov.matrix = cov.matrix))
  stat.dist [,1:6] <- sapply (matrix.stat.functions, MatrixStatisticsMap)
  stat.dist[,7] <- as.numeric (stat.dist[,5] > crit.value)
  stat.dist[,8] <- as.numeric (stat.dist[,6] > crit.value)
  stat.dist <- cbind (stat.dist, null.dist)
  colnames (stat.dist) <- c('respondability',
                            'evolvability',
                            'conditional.evolvability',
                            'autonomy',
                            'flexibility',
                            'constraints',
                            'flex.n',
                            'const.n',
                            'null.dist')
  stat.mean <- colMeans (stat.dist)
  stat.mean[7:8] <- stat.mean[7:8] * iterations
  pc1 <- eigen (cov.matrix)$vectors[,1]
  MatrixStatisticsMapPc1 <- function (CurrentFunc) return (CurrentFunc (beta = pc1, cov.matrix = cov.matrix))
  maximum <- sapply (matrix.stat.functions, MatrixStatisticsMapPc1)
  integration <- c (CalcR2 (cov.matrix), Pc1Percent (cov.matrix))
  names (integration) <- c ('MeanSquaredCorrelation', 'pc1%')
  stat.mean <- c (integration, stat.mean)
  return (list ('dist' = stat.dist, 'mean' = stat.mean, 'max.val' = maximum))
}
