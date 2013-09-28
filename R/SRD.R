SRD <- function (cov.matrix.1, cov.matrix.2, iterations = 1000)
  # Calculates the selection response decomposition comparison between covariance matrices
  #
  # Args:
  #     cov.matrix.(1,2): covariance matrices being compared
  #     iterations: number of RandomSkewers random vectors
  # Return:
  #     SRD scores for each trait and significance using mean and sd of SRD scores
{
  size <- dim (cov.matrix.1)[1]
  r2s <- array (0, c(size,iterations))
  beta <- apply (array (rnorm (size*iterations, mean = 0, sd = 1),c(size,iterations)),2, Normalize)
  for (I in 1:iterations){
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

PlotSRD <- function (output, matrix.label = "")
  # Plots the output of the SRD function in standard format
  #
  # Args:
  # output: the output from the SRD funtion
  # matrix.label: string with the names of the matrices that were compared in the SRD function
  # Return:
  # pretty plot
{
  layout (array (c(1,1,2,2),c(2,2)))
  par (mar = c(4.0, 4.0, 8, 0.4))
  mean.r2 <- output$out[,1]
  low.r2 <- output$out[,2]
  up.r2 <- output$out[,3]
  c.mean.r2 <- output$out[,5]
  c.sd.r2 <- output$out[,6]
  if (is.null (rownames (output$out)))
    dists <- 1:length (mean.r2)
  else
    dists <- rownames (output$out)
  ### plot scores
  ### b l t r
  pc.pch <- output$model$code + 20
  plot (mean.r2, type = "p", lty = 2, pch = pc.pch,
        ylab = "", xlab = "", xaxt = "n", ylim = c(-1,1))
  for (i in 1:length (mean.r2))
  {
    abline (v = i, lty = 1, col = rgb (0.8,0.8,0.8))
  }
  arrows (x0 = 1:length (mean.r2),
          y0 = low.r2, y1 = up.r2,
          angle = 90, length = 0.05, code = 3)
  abline (h = mean (mean.r2), lty = 3)
  axis (3, 1:length (mean.r2), dists, las = 2, cex.axis = 1.3)
  mtext (side = 2, at = mean (mean.r2), text = round (mean (mean.r2), 2), las = 2, cex=1.8)
  ### plot av sd
  ### b l t r
  par (mar = c(4.0, 0.0, 8, 4.6))

  output$model$code[output$model$code == 1] <- 1.2
  output$model$code[output$model$code == 0] <- 1
  pc.cex <- output$model$code
  plot (c.sd.r2 ~ c.mean.r2, pch = pc.pch, xlab = "", ylab = "",
        yaxt = "n", main = matrix.label, cex.main= 3, cex.axis=1.3)
  abline (v = 0, lty = 2)
  abline (h = 0, lty = 2)
  text (c.mean.r2, c.sd.r2, labels = dists, pos = 4, cex = pc.cex)
  axis (4, las = 2)
}
