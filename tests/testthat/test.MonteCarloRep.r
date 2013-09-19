test_that("MonteCarloRep returns sensible results",
          {
            expect_that(MonteCarloRep(diag(rep(2, 10)), 10) <= 1, is_true())
          }
)
