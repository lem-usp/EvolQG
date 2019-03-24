test_that("MINT returns correct results",{
  #skip_on_cran()
  modules = matrix(c(rep(c(1, 0, 0), each = 5),
                     rep(c(0, 1, 0), each = 5),
                     rep(c(0, 0, 1), each = 5)), 15)
  cor.hypot = CreateHypotMatrix(modules)[[4]]
  hypot.mask = matrix(as.logical(cor.hypot), 15, 15)
  mod.cor = matrix(NA, 15, 15)
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(42)
  mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
  mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
  diag(mod.cor) = 1
  mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
  hypothetical.modules = cbind(modules, matrix(sample(c(1, 0), 4*15, replace=TRUE), 15, 4))
  colnames(hypothetical.modules) <- letters[1:7]
  results = MINT(mod.cor, hypothetical.modules)
  #expect_equal(results[[1]], expected)
  expect_true(all(results[[2]] == modules))
  
  random_var = runif(15, 1, 10)
  out = raply(50, function(x){
    mod.cov = cov(mvtnorm::rmvnorm(50, sigma = sqrt(outer(random_var, random_var)) * mod.cor))
    hypothetical.modules = cbind(modules, matrix(sample(c(1, 0), 4*15, replace=TRUE), 15, 4))
    results = MINT(mod.cov, hypothetical.modules)[[1]]
    return(results[1,1])})
  close_enough = sum(grepl("1_2_3|2_3|1_2|1_3", out, perl = TRUE))/50
  expect_true(close_enough > 0.95)
  
  xdata = mvtnorm::rmvnorm(200, sigma = sqrt(outer(random_var, random_var)) * mod.cor)
  out_jack = JackKnifeMINT(xdata, hypothetical.modules, 100)
  close_enough_jack = sum(grepl("a_b_c", out_jack, perl = TRUE))/100
  expect_equal(out_jack[1,1], "a_b_c")
  expect_true(out_jack[1,5] > 0.95)
  
  expect_error(JackKnifeMINT(xdata, hypothetical.modules[,1]))
  expect_error(JackKnifeMINT(cov(xdata), hypothetical.modules))
})
