test_that("KrzCor returns correct results",
          {
              cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
              cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
              ret.dim = dim(cov.matrix.1)[1]/2 - 1
              EigenVectors <- function (x) return (eigen(x)$vectors[,1:ret.dim])
              A <- EigenVectors (cov.matrix.1)
              B <- EigenVectors (cov.matrix.2)
              S <- t(A) %*% B %*% t(B) %*% A
              SL <- sum (eigen(S)$values) / ret.dim
              expect_that(KrzCor(cov.matrix.1, cov.matrix.2), equals(SL))
              expect_that(KrzCor(cov.matrix.1, cov.matrix.2, 10), equals(1))
              expect_that(KrzCor(cov.matrix.1, cov.matrix.2) <= 1, is_true())
              expect_that(KrzCor(cov.matrix.1, cov.matrix.2) > 0, is_true())
          }
)
