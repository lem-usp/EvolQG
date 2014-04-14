test_that("BootstrapRep returns sensible results",
          {
              ind.data <- iris[,1:4]
              expect_that(BootstrapRep(ind.data, "randomskewers", 50) <= 1 , is_true())
              expect_that(BootstrapRep(ind.data, "mantel", 50) <= 1 , is_true())
              expect_that(BootstrapRep(ind.data, "krz", 50) <= 1 , is_true())
              expect_that(BootstrapRep(ind.data, "krz", 50, T) <= 1 , is_true())
              expect_that(BootstrapRep(ind.data, "randomskewers", 50) >= 0.98 , is_true())
              expect_that(BootstrapRep(ind.data, "randomskewers", 50, T) >= 0.98 , is_true())
              expect_that(BootstrapRep(ind.data, "mantel", 50) >= 0.98 , is_true())
              expect_that(BootstrapRep(ind.data, "krz", 50) >= 0.98 , is_true())
              expect_that(BootstrapRep(ind.data, "krz", 50, T) >= 0.98 , is_true())
              expect_that(BootstrapRep(cov(ind.data), "randomskewers", 50), throws_error("input appears to be a matrix, use residuals."))
              expect_that(BootstrapRep(cor(ind.data), "randomskewers", 50), throws_error("input appears to be a matrix, use residuals."))
          }
)
