MultiKrzCor <-
function (matrix.list, repeat.vector = NULL)
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
  n.matrix <- length (matrix.list)
  matrix.names <- names (matrix.list)
  correlations <- array (0, c(n.matrix, n.matrix))
  for (i in 1:(n.matrix - 1)) {
    for (j in (i+1):n.matrix) {
      cat (i, ' ', j, '\n')
      comparing.now <- KrzCor (matrix.list [[i]],
                               matrix.list [[j]])
      correlations [i, j] <- comparing.now
      if (!is.null (repeat.vector))
        correlations [j, i] <- correlations [i, j] / sqrt (repeat.vector [i] * repeat.vector [j])
    }
  }
  if (!is.null (repeat.vector)) {
    diag (correlations) <- repeat.vector
  }
  rownames (correlations) <- matrix.names
  colnames (correlations) <- matrix.names
  return (correlations)
}
