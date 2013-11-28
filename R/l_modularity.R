#' L Modularity
#' 
#' Calculates the L-Modularity and a partition of traits
#' @export
#' @useDynLib Morphometrics
#' @examples
#' cor.matrix = RandomCorr(10)
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
  return(list(l_modularity, mod_hipotesis))
}