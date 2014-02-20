#' L Modularity
#' 
#' Calculates the L-Modularity and a partition of traits
#' @param cor.matrix correlation matrix
#' @return List with L-Modulariy value and trait partition
#' @export
#' @useDynLib Morphometrics
#' @importFrom Rcpp evalCpp
#' @examples
#' cor.matrix = RandomMatrix(10)
#' LModularity(cor.matrix)
LModularity <- function(cor.matrix){
  num_traits <- dim(cor.matrix)[1]
  s <- numeric(num_traits)
  l_modularity <- annealing(cor.matrix, s)
  partition <- as.numeric(factor(s))
  modules <- unique(partition)
  num_modules <- length(modules)
  mod_hipotesis <- array(0, c(num_traits, num_modules))
  for (mod in modules){
    mod_hipotesis[partition == mod, mod] = 1
  }
  output = list("L-Modularity" = l_modularity, "Modularity Hipotesis" = mod_hipotesis) 
  return(output)
}