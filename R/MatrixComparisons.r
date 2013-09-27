ComparisonMap <- function (matrix.list, MatrixCompFunc, repeat.vector = NULL, num.cores = 1){
  # Performs multiple comparisons between a set of covariance or
  # correlation matrices.
  #
  # Args:
  #  matrix.list: a list of covariance or correlation matrices
  #  MatrixCompFunc: function to use for comparison
  #  repeat.vector: vector of matrix repeatabilities
  #
  # Return:
  #  a list with two matrices containing $\Gamma$-values or average random
  #  skewers correlation and probabilities according to permutation test.
  #  if repeat.vector was also passed, values below the diagonal on the correlation matrix
  #  will contain corrected correlation values.
  library(plyr)
  library(reshape2)
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else{
    parallel = FALSE
  }
  n.matrix <- length(matrix.list)
  if(is.null(names(matrix.list))) {names(matrix.list) <- 1:n.matrix}
  matrix.names <- names (matrix.list)
  CompareToN <- function(n) ldply(matrix.list[(n+1):n.matrix],
                                  function(x) {MatrixCompFunc(x, matrix.list[[n]])[1:2]},
                                  .parallel = parallel)
  comparisons <- adply(1:(n.matrix-1), 1,  CompareToN, .parallel = parallel)
  corrs <- acast(comparisons[-4], X1~.id)[,matrix.names[-1]]
  probs <- acast(comparisons[-3], X1~.id)[,matrix.names[-1]]
  probabilities <- array (0, c(n.matrix, n.matrix))
  correlations <- probabilities
  probabilities[upper.tri(probabilities)] <- probs[upper.tri(probs, diag=T)]
  correlations[upper.tri(correlations)] <- corrs[upper.tri(probs, diag=T)]
  if (!is.null (repeat.vector)) {
    repeat.matrix <- sqrt(outer(repeat.vector, repeat.vector))
    correlations[lower.tri(correlations)] <- t(correlations/repeat.matrix)[lower.tri(correlations)]
    diag (correlations) <- repeat.vector
  }
  rownames (correlations) <- matrix.names
  colnames (correlations) <- matrix.names
  dimnames (probabilities) <- dimnames (correlations)
  output <- list ('correlations' = correlations, 'probabilities' = probabilities)
  return (output)
}

SingleComparisonMap  <- function(matrix.list, y.mat, MatrixCompFunc, num.cores){
  library(plyr)
  library(reshape2)
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else
    parallel = FALSE
  if(!all.equal(dim(matrix.list[[1]]),dim(y.mat)))
    stop("Matrices on list and single matrice dimension do not match")
  else
    out <- ldply(matrix.list,
                 function(x) {MatrixCompFunc(x, y.mat)},
                 .parallel = parallel)
  return(out)
}

RandomSkewers <- function(cov.matrix, ...) UseMethod("RandomSkewers")

RandomSkewers.default <- function (cov.matrix.1, cov.matrix.2, iterations = 1000)
  # Calculates covariance matrix correlation via random skewers
  # Args:
  #     cov.matrix.(1,2): Two covariance matrices to be compared
  #     iterations: Number of generated random skewers
  # Return:
  #     List with mean value of correlation, p value and standard deviation
{
  traits <- dim (cov.matrix.1) [1]
  base.vector <- Normalize(rnorm(traits))
  random.vectors <- array (rnorm (iterations * traits, mean = 0, sd = 1), c(traits, iterations))
  random.vectors <- apply (random.vectors, 2, Normalize)
  dist <- base.vector %*% random.vectors
  dz1 <- apply (cov.matrix.1 %*% random.vectors, 2, Normalize)
  dz2 <- apply (cov.matrix.2 %*% random.vectors, 2, Normalize)
  real <- apply (dz1 * dz2, 2, sum)
  ac <- mean (real)
  stdev <- sd (real)
  prob <- sum (ac < dist) / iterations
  output <- c(ac, prob, stdev)
  names(output) <- c("AC","P","SD")
  return(output)
}

RandomSkewers.list <- function (matrix.list, cov.matrix.2 = NULL, repeat.vector = NULL, iterations = 1000, num.cores = 1)
{
  if (is.null (cov.matrix.2)) {
    out <- ComparisonMap(matrix.list,
                         function(x, y) RandomSkewers.default(x, y, iterations),
                         repeat.vector = repeat.vector,
                         num.cores = num.cores)
  }
  else{
    out <- SingleComparisonMap(matrix.list, cov.matrix.2,
                               function(x, y) RandomSkewers.default(x, y, iterations),
                               num.cores = num.cores)
  }
  return(out)
}

MantelCor <- function (cor.matrix, ...) UseMethod("MantelCor")

MantelCor.default <- function (cor.matrix.1, cor.matrix.2, iterations = 100, mod = FALSE)
  # Calculates matrix correlation with confidence intervals using mantel permutations
  #
  # Args:
  #     cor.matrix.(1,2): correlation matrices being compared
  #     iterations: number of permutations
  #     mod: for when testing binary modularity hipotesis
  # Return:
  #     matrix pearson correelation and significance.
  #     if mod==TRUE also returns average within, between and average ratio correlations
{
  library(vegan)
  mantel.out <- mantel(cor.matrix.1, cor.matrix.2, permutations = iterations)
  correlation <- mantel.out$statistic
  prob <- mantel.out$signif
  if (mod == TRUE){
    index <- cor.matrix.2[lower.tri(cor.matrix.2)]
    avg.plus <- mean (cor.matrix.1 [lower.tri(cor.matrix.1)] [index != 0])
    avg.minus <- mean (cor.matrix.1 [lower.tri(cor.matrix.1)] [index == 0])
    avg.ratio <- avg.plus / avg.minus
    output <- c(correlation, prob, avg.plus, avg.minus, avg.ratio)
    names(output) <- c("Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio")
  }
  else{
    output <- c(correlation, prob)
    names(output) <- c("Rsquared", "Probability")
  }
  return (output)
}

MantelCor.list <- function (matrix.list, cor.matrix.2 = NULL,
                            repeat.vector = NULL, iterations = 100, mod = FALSE, num.cores = 1)
{
  if (is.null (cor.matrix.2)) {
    out <- ComparisonMap(matrix.list,
                         function(x, cor.matrix.2) MantelCor.default(x, cor.matrix.2, iterations),
                         repeat.vector = repeat.vector,
                         num.cores = num.cores)
  }
  else{
    out <- SingleComparisonMap(matrix.list, cor.matrix.2,
                               function(x, y) MantelCor.default(y,
                                                                         x,
                                                                         iterations, mod = mod),
                               num.cores = num.cores)
  }
  return(out)
}

KrzCor <- function (cov.matrix, ...) UseMethod("KrzCor")

KrzCor.default <- function (cov.matrix.1, cov.matrix.2, ret.dim = NULL)
  # Calculates the Krzanowski correlation between matrices
  #
  # Args:
  #     cov.matrix.(1,2): covariance being compared
  #     ret.dim: number of retained dimensions in the comparison,
  #              default for nxn matrix is n/2-1
  # Return:
  #     Kzranowski correlation
{
  if (is.null(ret.dim))
    ret.dim = round(dim(cov.matrix.1)[1]/2 - 1)
  EigenVectors <- function (x) return (eigen(x)$vectors[,1:ret.dim])
  A <- EigenVectors (cov.matrix.1)
  B <- EigenVectors (cov.matrix.2)
  S <- t(A) %*% B %*% t(B) %*% A
  SL <- sum (eigen(S)$values) / ret.dim
  return (SL)
}

KrzCor.list <- function (matrix.list, cov.matrix.2 = NULL,
                         repeat.vector = NULL, ret.dim = NULL,
                         num.cores = 1)
  # Performs multiple comparisons between a set of covariance or
  # correlation matrices using Kzranowski Correlation.
  #
  # Args:
  #  matrix.list: a list of covariance or correlation matrices
  #  repeat.vector: vector of matrix repeatabilities
  #
  # Return:
  #  Matrices containing Krzanowski correaltion between matrices.
  #  if repeat.vector was also passed, values below the diagonal on the correlation matrix
  #  will contain corrected correlation values.
{
  if (is.null (cov.matrix.2)) {
    out <- ComparisonMap(matrix.list,
                         function(x, y) return(c(KrzCor.default(x, y, ret.dim), NA)),
                         repeat.vector = repeat.vector,
                         num.cores = num.cores)
    out <- out[[1]]
  }
  else{
    out <- SingleComparisonMap(matrix.list, cov.matrix.2,
                         function(x, y) return(c(KrzCor.default(x, y, ret.dim), NA)),
                               num.cores = num.cores)
    out <- out[,-length(out)]
  }
  return(out)
}
