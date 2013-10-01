AlphaRep <- function (cor.matrix, tam)
  # Calculates the matrix repeatability using the equation in Cheverud 1996
  # Quantitative genetic analysis of cranial morphology in the cotton-top
  # (Saguinus oedipus) and saddle-back (S. fuscicollis) tamarins. Journal of Evolutionary Biology 9, 5-42.
  #
  # Args:
  #     cor.matrix: correlation matrix
  #     tam: sample size
  # Return:
  #     matrix repeatability
{
  vec <- cor.matrix[lower.tri(cor.matrix)]
  var.erro <- (1 - mean(vec)^2)/(tam-2)
  var.vec <- var(vec)
  return((var.vec - var.erro)/var.vec)
}
