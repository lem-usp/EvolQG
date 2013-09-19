test_that("MonteCarloR2 returns sensible results",
          {
            expect_that(length(MonteCarloR2(diag(rep(1, 10)), 10)), equals(1000))
            expect_that(length(MonteCarloR2(diag(rep(1, 10)), 10, 10)), equals(10))
            results <- MonteCarloR2(diag(rep(1, 10)), 10)
            corrs = sapply(results, function(x) isTRUE(x < 1 & x > 0))
            expect_that(sum(corrs), equals(1000))
          }
)
