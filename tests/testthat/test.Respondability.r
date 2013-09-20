test_that("Respondability returns correct result",
         {
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           repond <- Norm(cov.matrix %*% beta)
           expect_that(Respondability(beta, cov.matrix), equals(repond))
         }
)

