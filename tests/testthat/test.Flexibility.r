test_that("Flexibility returns correct result",
         {
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           flex <- t (beta) %*% cov.matrix %*% beta / Norm (cov.matrix %*% beta)
           expect_that(Flexibility(beta, cov.matrix), equals(flex))
           expect_that(Flexibility(beta, cov.matrix) <=  1, is_true())
           expect_that(Flexibility(beta, cov.matrix) >= -1, is_true())
         }
)

