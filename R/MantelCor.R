MantelCor <-
function (cor.matrix.1, cor.matrix.2, nit = 1000, mod = FALSE)
  # Calculates matrix correlation with confidence intervals using mantel permutations
  #
  # Args:
  #     cor.matrix.(1,2): correlation matrices being compared
  #     nit: number of permutations
  #     mod: for when testing binary modularity hipotesis
  # Return:
  #     matrix pearson correelation and significance.
  #     if mod==TRUE also returns average within, between and average ratio correlations
{
  if(!require(vegan)) install.packages("vegan")
  library(vegan)
  mantel.out <- mantel(cor.matrix.1, cor.matrix.2)
  correlation <- mantel.out$statistic
  prob <- mantel.out$signif
  if (mod == TRUE){
    avg.plus <- mean (cor.matrix.1 [lower.tri(cor.matrix.1)] [index != 0])
    avg.minus <- mean (cor.matrix.1 [lower.tri(cor.matrix.1)] [index == 0])
    avg.ratio <- avg.plus / avg.minus
    output <- c(correlation, prob, avg.plus, avg.minus, avg.ratio)
    names(output) <- c("R²", "Probability", "AVG+", "AVG-", "AVG Ratio")
  }
  else{
    output <- c(correlation, prob)
    names(output) <- c("R²", "Probability")
  }
  return (output)
}
