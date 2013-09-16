PlotSRD <-
function (output, matrix.label = "")
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
