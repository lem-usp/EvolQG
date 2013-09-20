MultiRsMantel <-
function (matrix.list, MatrixCompFunc = RandomSkewers, repeat.vector = NULL, iterations = 1000)
  # Performs multiple comparisons between a set of covariance or
  # correlation matrices.
  #
  # Args:
  #  matrix.list: a list of covariance or correlation matrices
  #  MatrixCompFunc: function to use for comparison
  #  repeat.vector: vector of matrix repeatabilities
  #  iterations: number of RandomSkewers or matrix permutations passed to MatrixCompFunc
  #
  # Return:
  #  a list with two matrices containing $\Gamma$-values or average random
  #  skewers correlation and probabilities according to permutation test.
  #  if repeat.vector was also passed, values below the diagonal on the correlation matrix
  #  will contain corrected correlation values.
{
  n.matrix <- length (matrix.list)
  matrix.names <- names (matrix.list)
  probabilities <- array (0, c(n.matrix, n.matrix))
  correlations <- probabilities
  for (i in 1:(n.matrix - 1)) {
    for (j in (i+1):n.matrix) {
      comparing.now <- MatrixCompFunc (matrix.list [[i]],
                                       matrix.list [[j]],
                                       iterations)
      correlations [i, j] <- comparing.now [1]
      probabilities [i, j] <- comparing.now [2]
      if (!is.null (repeat.vector))
        correlations [j, i] <- correlations [i, j] / sqrt (repeat.vector [i] * repeat.vector [j])
    }
  }
  if (!is.null (repeat.vector)) {
    diag (correlations) <- repeat.vector
  }
  rownames (correlations) <- matrix.names
  colnames (correlations) <- matrix.names
  dimnames (probabilities) <- dimnames (correlations)
  output <- list ('correlations' = correlations, 'probabilities' = probabilities)
  return (output)
}
