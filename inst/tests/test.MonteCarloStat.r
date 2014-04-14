test_that("MonteCarloR2 returns sensible results",
          {
            expect_that(length(MonteCarloR2(diag(rep(2, 10)), 10)), equals(1000))
            expect_that(length(MonteCarloR2(diag(rep(2, 10)), 10, 10)), equals(10))
            results <- MonteCarloR2(diag(rep(2, 10)), 10)
            corrs = sapply(results, function(x) isTRUE(x < 1 & x > 0))
            expect_that(sum(corrs), equals(1000))
            expect_that(MonteCarloR2(RandomMatrix(10), 10), gives_warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling."))
          }
)
test_that("MonteCarloRep returns sensible results",
          {
            cov.matrix <- RandomMatrix(10, 1, 1, 10)
            expect_that(MonteCarloRepRandomSkewers(RandomMatrix(10), 30), gives_warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling."))
            expect_that(MonteCarloRepRandomSkewers(cov.matrix, 30) <= 1, is_true())
            expect_that(MonteCarloRepMantelCor(cov.matrix, 30) <= 1, is_true())
            expect_that(MonteCarloRepKrzCor(cov.matrix, 30) <= 1, is_true())
            expect_that(MonteCarloRepKrzCor(cov.matrix, 30, T) <= 1, is_true())
          }
)
