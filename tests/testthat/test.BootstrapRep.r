test_that("BootstrapRep returns sensible results",
          {
              ind.data <- iris[,1:4]
              expect_that(BootstrapRepRandomSkewers(ind.data, 50) <= 1 , is_true())
              expect_that(BootstrapRepMantelCor(ind.data, 50) <= 1 , is_true())
              expect_that(BootstrapRepKrzCor(ind.data, 50) <= 1 , is_true())
              expect_that(BootstrapRepKrzCor(ind.data, 50, T) <= 1 , is_true())
              expect_that(BootstrapRepRandomSkewers(ind.data, 50) >= 0.98 , is_true())
              expect_that(BootstrapRepRandomSkewers(ind.data, 50, T) >= 0.98 , is_true())
              expect_that(BootstrapRepMantelCor(ind.data, 50) >= 0.98 , is_true())
              expect_that(BootstrapRepKrzCor(ind.data, 50) >= 0.98 , is_true())
              expect_that(BootstrapRepKrzCor(ind.data, 50, T) >= 0.98 , is_true())
              expect_that(BootstrapRepRandomSkewers(cov(ind.data), 50), throws_error("input appears to be a matrix, use residuals."))
              expect_that(BootstrapRepRandomSkewers(cor(ind.data), 50), throws_error("input appears to be a matrix, use residuals."))
          }
)
