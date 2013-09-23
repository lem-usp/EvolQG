TestModularity <- function (cor.matrix, modularity.hipot)
  # Tests modularity hipotesis using cor.matrix matrix and trait groupings
  #
  # Args:
  #   cor.matrix: Correlation matrix
  #   modularity.hipot: Matrix of hipotesis.
  #                     Each line represents a trait and each column a module.
  #                     if modularity.hipot[i,j] == 1, trait i is in module j.
  # Return:
  #   list with mantel correlation of cor.matrix with binary hipotesis matrices
{
  no.hip <- dim (modularity.hipot) [2] + 1
  m.hip.array <- CreateHipotMatrix(modularity.hipot)
  output <- array (0, c(5, no.hip))
  for (N in 1:no.hip){
    tmp <- MantelCor (cor.matrix,
                      m.hip.array[,,N],
                      mod = TRUE)
    output[,N] <- tmp
  }
  if(is.null(colnames(modularity.hipot))) colnames(modularity.hipot) <- 1:(no.hip-1)
  dimnames(output) <- list (names (tmp), c (colnames (modularity.hipot),"Full Integration"))
  return (output)
}

CreateHipotMatrix <- function(modularity.hipot)
{
  no.hip <- dim (modularity.hipot) [2]
  traits <- dim (modularity.hipot) [1]
  m.hip.array <- array (0, c(traits, traits, no.hip + 1))
  for (N in 1:no.hip){
    for (L in 1:traits){
      for (M in 1:traits){
        m.hip.array[L,M,N] <- ifelse (modularity.hipot[L,N] &
                                      modularity.hipot[M,N],
                                      1,
                                      0)
      }
    }
  }
  m.hip.array[,,no.hip+1] <- as.integer (as.logical (apply (m.hip.array, c(1,2), sum)))
  return(m.hip.array)
}
