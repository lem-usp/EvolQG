CalculateMatrix <-
function(linear.m){
  # Calculates the covariance matrix using the residuals of a linear model
  #
  # Args:
  #   linear.m: a linear model created with the lm() function
  # Return:
  #   cov.matrix: the covariance matrix of the residuals
  cov.matrix = var(linear.m$residuals)*((dim(linear.m$residuals)[1]-1)/linear.m$df.residual)
  return (cov.matrix)
}
