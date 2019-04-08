test_that("LModularity resturn resonable results",{
  skip_on_cran()
  modules = matrix(c(rep(c(1, 0, 0), each = 5),
                     rep(c(0, 1, 0), each = 5),
                     rep(c(0, 0, 1), each = 5)), 15)
  cor.hypot = CreateHypotMatrix(modules)[[4]]
  hypot.mask = matrix(as.logical(cor.hypot), 15, 15)
  mod.cor = matrix(NA, 15, 15)
  mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
  mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
  diag(mod.cor) = 1
  mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
  expect_equal(modules, LModularity(mod.cor)[[2]])
  expect_warning(squared <- LModularity(-mod.cor))
  expect_equivalent(squared[[2]], modules)
})
