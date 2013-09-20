test_that("Autonomy returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              auto <- ((t (beta) %*% solve (cov.matrix) %*% beta)^(-1)) / (t (beta) %*% cov.matrix %*% beta)
              expect_that(Autonomy(beta, cov.matrix), equals(auto))
          }
)
