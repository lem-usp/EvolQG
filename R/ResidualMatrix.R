ResidualMatrix <-
function (model)
  {
    ## Calculates residual matrix from a estimated linear model
    ## Args:
    ##  model: linear model previously estimated
    ## Value:
    ##  cov.mat: residual covariance matrix
    res <- residuals (model)
    res.df <- model $ df.residual
    cov.mat <- t (res) %*% res / res.df
    return (cov.mat)
  }
