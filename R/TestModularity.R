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
  m.hip.list <- CreateHipotMatrix(modularity.hipot)
  if(is.null(colnames(modularity.hipot))) colnames(modularity.hipot) <- 1:dim (modularity.hipot) [2]
  names(m.hip.list) <- c(colnames (modularity.hipot),"Full Integration")
  output <- MantelCor (m.hip.list, cor.matrix, mod = TRUE)
  return (output)
}

CreateHipotMatrix <- function(modularity.hipot)
{
  num.hip <- dim (modularity.hipot) [2]
  num.traits <- dim (modularity.hipot) [1]
  m.hip.list <- alply(modularity.hipot, 2, function(x) outer(x, x))
  m.hip.list[[num.hip+1]] <- matrix(as.integer (as.logical (Reduce ("+", m.hip.list[1:num.hip]))),
                                    num.traits, num.traits, byrow=T)
  return(m.hip.list[1:(num.hip+1)])
}
