test_that("Constraints returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta = rnorm(4)
              const = abs (t (Normalize (eigen (cov.matrix)$vectors[,1])) %*% Normalize (cov.matrix %*% beta))
              expect_that(Constraints(beta, cov.matrix), equals(const))
          }
)
