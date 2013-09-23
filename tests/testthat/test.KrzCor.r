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

test_that("KrzCor returns correct results",
          {
            cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.3 <- cov(matrix(rnorm(30*10), 30, 10))
            mat.list <- list(cov.matrix.1, cov.matrix.2, cov.matrix.3)
            results <- KrzCor(mat.list)
            results.2 <- KrzCor(mat.list, repeat.vector = c(0.8, 0.9, 0.85))
            expect_that(sum(is.na(results)), equals(0))
            expect_that(sum(is.na(results.2)), equals(0))
            expect_that(dim(results), equals(c(length(mat.list),length(mat.list))))
            upper  <- results[upper.tri(results)]
            upper.bool = sapply(upper, function(x) isTRUE(x > 0 & x < 1))
            expect_that(sum(upper.bool), equals(length(upper.bool)))
            lower  <- results[lower.tri(results, diag = T)]
            lower.bool = sapply(lower, function(x) isTRUE(x == 0))
            expect_that(sum(lower.bool), equals(length(lower.bool)))
            upper.2  <- results.2[upper.tri(results.2)]
            lower.2  <- results.2[lower.tri(results.2)]
            expect_that(upper.2, equals(upper))
            expect_that(diag(results.2), is_equivalent_to(c(0.8, 0.9, 0.85)))
            expect_that(sum(lower.2 > upper.2), equals(length(mat.list)))
          }
)
