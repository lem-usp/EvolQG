SRD <-
function (cov.matrix.1, cov.matrix.2, nsk = 1000)
  # Calculates the selection response decomposition comparison between covariance matrices
  #
  # Args:
  #     cov.matrix.(1,2): covariance matrices being compared
  #     nsk: number of RandomSkewers random vectors
  # Return:
  #     SRD scores for each trait and significance using mean and sd of SRD scores
{
  size <- dim (cov.matrix.1)[1]
  r2s <- array (0, c(size,nsk))
  beta <- apply (array (rnorm (size*nsk, mean = 0, sd = 1),c(size,nsk)),2, Normalize)
  for (I in 1:nsk){
    beta.matrix <- diag (beta[,I])
    dz1 <- apply (cov.matrix.1 %*% beta.matrix, 1, Normalize)
    dz2 <- apply (cov.matrix.2 %*% beta.matrix, 1, Normalize)
    r2s[,I] <- colSums (dz1 * dz2)
  }
  # results
  mean.r2 <- apply (r2s, 1, mean)
  sample.conf <- function (x, lower = TRUE){
    ox <- x[order(x)]
    lox <- length (ox)
    if (lower)
      crit <- round (0.025 * lox)
    else
      crit <- round (0.975 * lox)
    return (ox[crit])
  }
  low.r2 <- apply (r2s, 1, sample.conf, lower = TRUE)
  up.r2 <- apply (r2s, 1, sample.conf, lower = FALSE)
  sd.r2 <- apply (r2s,1,sd)
  cmean.r2 <- scale (mean.r2, scale = FALSE)
  csd.r2 <- scale (sd.r2, scale = FALSE)
  cent <- cbind (cmean.r2,csd.r2)
  pca.cent <- princomp (cent, cor = TRUE, scores = TRUE)
  pc1 <- pca.cent$scores[,1]
  if (pca.cent$loadings[1,1] < 0)
    pc1 <- - pc1
  pc1 <- pc1 / sd (pc1)
  pc1.quant <- quantile (pc1,
                         probs = c (1,5,10,20,25,30,40,50,60,70,75,80,90,95,99)/100,
                         names = FALSE)
  pc1.int <- - 1.96 / sqrt (length (pc1))
  pc1.sig <- ifelse (pc1 < pc1.int, 1, 0)
  model <- list ("quantiles" = pc1.quant,
                 "interval"  = pc1.int,
                 "code"      = pc1.sig)
  output <- cbind (mean.r2, low.r2, up.r2, sd.r2, cmean.r2, csd.r2)
  colnames (output) <- c("ARC","IC-","IC+","SD","CMEAN","CSD")
  rownames (output) <- rownames (cov.matrix.1)
  return (list ("out"    = output,
                "pc1"    = pc1,
                "model"  = model,
                "cormat" = cor (t(r2s))))
}
