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
            expect_that(MonteCarloRep(RandomMatrix(10), RandomSkewers, 30), gives_warning("Matrix appears to be a correlation matrix! Only covariance matrices should be used in parametric resampling."))
            expect_that(MonteCarloRep(cov.matrix, RandomSkewers, 30) <= 1, is_true())
            expect_that(MonteCarloRep(cov.matrix, MantelCor, 30) <= 1, is_true())
            expect_that(MonteCarloRep(cov.matrix, KrzCor, 30) <= 1, is_true())
            expect_that(MonteCarloRep(cov.matrix, KrzCor, 30, T) <= 1, is_true())
          }
)
test_that("MonteCarloStat throws error",
{
  expect_that(MonteCarloStat(array(1:100, c(10, 10)), 10, 10, RandomSkewers, cov), 
              throws_error("covariance matrix must be symmetric."))
})
