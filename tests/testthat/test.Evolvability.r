test_that("Evolvability returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- rnorm(4)
              evol <- t (beta) %*% cov.matrix %*% beta
              expect_that(Evolvability(beta, cov.matrix), equals(evol))
          }
)
