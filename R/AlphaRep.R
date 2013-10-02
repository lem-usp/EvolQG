AlphaRep <- function (cor.matrix, sample.size){
  vec <- cor.matrix[lower.tri(cor.matrix)]
  var.erro <- (1 - mean(vec)^2)/(sample.size-2)
  var.vec <- var(vec)
  return((var.vec - var.erro)/var.vec)
}
