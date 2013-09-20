test_that("Pc1Percent returns correct results",
          {
            cov.matrix = cov(matrix(rnorm(30*10), 30, 10))
            pc1 <- eigen(cov.matrix)$values[1]/sum(diag(cov.matrix))
            expect_that(Pc1Percent(cov.matrix), equals(pc1))
          }
)
