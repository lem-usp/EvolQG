CalcR2CvCorrected  <- function(ind.data, ...) UseMethod("CalcR2CvCorrected")

CalcR2CvCorrected.default <- function (ind.data, iterations = 1000, num.cores = 1, ...) {
  cv <- function (x) return (sd(x)/mean(x))
  Stats = function(x) {
    cov.matrix = var(x)
    cor.matrix = cov2cor(cov.matrix)
    return(c(CalcR2(cor.matrix), cv(eigen(cov.matrix)$values), mean (apply (x, 2, cv))))
  }
  it.stats <- BootstrapRep(ind.data, iterations,
                           ComparisonFunc = function(x, y) y,
                           StatFunc = Stats,
                           num.cores = num.cores)
  colnames(it.stats) <- c("r2", "eVals_cv", "mean_cv")
  return (it.stats)
}

CalcR2CvCorrected.lm <- function (ind.data, iterations = 1000, ...) {
  cv <- function (x) return (sd(x)/mean(x))
  orv <- model$model[[1]]
  fac <- model$model[-1]
  df <- model$df.residual
  res <- residuals (model)
  size <- dim (res)
  straps <- array (0, c(size[1], iterations))
  r2 <- mcv <- evar <- c()
  i <- 1
  while (i <= iterations)
  {
    straps[,i] <- sample ((size[1]), replace = TRUE)
    corr <- cov2cor (var (res[straps[,i],]) * (size[1] - 1)/df)
    r2[i] <- mean (corr[lower.tri(corr)]^2)
    evar[i] <- cv (eigen(var (res[straps[,i],] * (size[1] - 1)/df))$values)
    orv.s <- orv[straps[,i],]
    fac.s <- fac[straps[,i],]
    tmp1 <- table (fac.s)
    tmp2 <- by (orv.s,fac.s,cv)
    tmp3 <- unlist (lapply (tmp2, mean))
    mcv[i] <- sum ((tmp3 * tmp1)/sum (tmp1))
    if (!is.na(mcv[i]))
      i <- i + 1
  }
  corr.or <- cov2cor (residual.matrix (model))
  r2.or <- mean (corr.or[lower.tri(corr.or)]^2)
  evar.or <- cv (eigen(residual.matrix (model))$values)
  tmp1 <- table (fac)
  tmp2 <- by (orv,fac,cv)
  tmp3 <- unlist (lapply (tmp2, mean))
  mcv.or <- sum ((tmp3 * tmp1)/sum (tmp1))
  return (list("sims" = cbind (r2, evar, mcv),
               "originals" = c(r2.or,evar.or,mcv.or),
               "CORmat" = corr.or,
               "perms" = straps))
}
