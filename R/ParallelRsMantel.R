ParallelRsMantel <-
function (matrix.list, MatrixCompFunc = RandomSkewers,
                              repeat.vector = NULL, iterations = 1000, n.cores = 2)
  ## Performs multiple comparisons between a set of covariance or
  ## correlation matrices.
  ##
  ## Args:
  ##  matrix.list: a list of covariance or correlation matrices
  ##  MatrixCompFunc: function to use for comparison
  ##  repeat.vector: vector of matrix repeatabilities
  ##  iterations: number of RandomSkewers or
  ##  matrix permutations passed to MatrixCompFunc
  ##
  ## Return:
  ##  a list with two matrices containing $\Gamma$-values or average random
  ##  skewers correlation and probabilities according to permutation test.
  ##  if repeat.vector was also passed, values below the diagonal
  ## on the correlation matrix
  ##  will contain corrected correlation values.
{
  require (multicore)
  n.matrix <- length (matrix.list)
  matrix.names <- names (matrix.list)
  correlations <- array (0, c(n.matrix, n.matrix))
  rownames (correlations) <- matrix.names
  colnames (correlations) <- matrix.names
  probabilities <- correlations
  index <- which (upper.tri (diag (n.matrix)), arr.ind = TRUE)
  singleComp = function (k)
    {
      i <- index [k,1]
      j <- index [k,2]
      cat (i, ' ', j, '\n')
      comparing.now <- MatrixCompFunc (matrix.list [[i]],
                                       matrix.list [[j]],
                                       iterations)
      return (comparing.now)
    }
  corr.list = mclapply (1:nrow (index), singleComp, mc.cores = n.cores)
  for (k in 1:length (corr.list))
    {
      i <- index [k,1]
      j <- index [k,2]
      correlations [i, j] <- corr.list [[k]] [1]
      probabilities [i, j] <- corr.list [[k]] [2]
      if (!is.null (repeat.vector))
        correlations [j, i] <- correlations [i, j] /
          sqrt (repeat.vector [i] * repeat.vector [j])
    }
  if (!is.null (repeat.vector)) {
    diag (correlations) <- repeat.vector
  }
  output <- list ('correlations' = correlations, 'probabilities' = probabilities)
  return (output)
}
