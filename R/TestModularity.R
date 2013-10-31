#' Test modularity hipotesis
#'
#' Tests modularity hipotesis using cor.matrix matrix and trait groupings
#' @param cor.matrix Correlation matrix
#' @param modularity.hipot modularity.hipot: Matrix of hipotesis.
#' Each line represents a trait and each column a module.
#' if modularity.hipot[i,j] == 1, trait i is in module j.
#' @param iterations Number of iterations, to be passed to MantelCor
#' @return Returns mantel correlation and associated probability for each modularity hipotesis, along with AVG+, AVG-, AVG Ratio for each module.
#' A total hipotesis combining all hipotesis is also tested.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{MantelCor}}
#' @export
#' @rdname TestModularity
#' @examples
#' cor.matrix <- RandomMatrix(10)
#' rand.hipots <- matrix(sample(c(1, 0), 30, replace=TRUE), 10, 3)
#' mod.test <- TestModularity(cor.matrix, rand.hipots)
#' @keywords mantel
#' @keywords modularity
TestModularity <- function (cor.matrix, modularity.hipot, iterations = 100) {
  m.hip.list <- CreateHipotMatrix(modularity.hipot)
  if(is.null(colnames(modularity.hipot))) colnames(modularity.hipot) <- 1:dim (modularity.hipot) [2]
  names(m.hip.list) <- c(colnames (modularity.hipot),"Full Integration")
  output <- MantelCor (m.hip.list, cor.matrix, iterations = iterations, mod = TRUE)
  return (output)
}

#' @export
#' @rdname TestModularity
CreateHipotMatrix <- function(modularity.hipot) {
  num.hip <- dim (modularity.hipot) [2]
  num.traits <- dim (modularity.hipot) [1]
  m.hip.list <- alply(modularity.hipot, 2, function(x) outer(x, x))
  m.hip.list[[num.hip+1]] <- matrix(as.integer (as.logical (Reduce ("+", m.hip.list[1:num.hip]))),
                                    num.traits, num.traits, byrow=T)
  return(m.hip.list[1:(num.hip+1)])
}
