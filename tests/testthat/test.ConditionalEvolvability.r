test_that("ConditionalEvolvability returns correct restults",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- rnorm(4)
              cond.evol <- (t (beta) %*% solve (cov.matrix) %*% beta)^(-1)
              expect_that(ConditionalEvolvability(beta, cov.matrix), equals(cond.evol))
          }
)
