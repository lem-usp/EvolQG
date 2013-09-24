test_that("MonteCarloRep returns sensible results",
          {
            cov.matrix <- RandomMatrix(10, 1, 1, 10)
            expect_that(MonteCarloRepRandomSkewers(cov.matrix, 30) <= 1, is_true())
            expect_that(MonteCarloRepMantelCor(cov.matrix, 30) <= 1, is_true())
            expect_that(MonteCarloRepKrzCor(cov.matrix, 30) <= 1, is_true())
            expect_that(MonteCarloRepKrzCor(cov.matrix, 30, T) <= 1, is_true())
          }
)
