CalcRepeatability <-
function (ID, ind.data)
  # Calculates Repeatabilities acording to:
  #    Lessels, C. M., & Boag, P. T. (1987).
  #    Unrepeatable repeatabilities: a common mistake.
  #    The Auk, 2(January), 116–121.
  # Args:
  #     ID: indentity of individuals
  #     ind.data: individual measurments
  # Return:
  #     vector of repeatabilities
{
  models.list <- apply (ind.data, 2, function (vec){return (lm (vec ~ ID))})
  models.list <- lapply (models.list, anova)
  rep.itself <- function (lm.model){
    msq <- lm.model$'Mean Sq' ## 1 entre, 2 dentro
    s2a <- (msq[1] - msq[2])/2
    out <- s2a / (s2a + msq[2])
    return (out)
  }
  out <- sapply (models.list, rep.itself)
  names (out) <- colnames (ind.data)
  return (out)
}
