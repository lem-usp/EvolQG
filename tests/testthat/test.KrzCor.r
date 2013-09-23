test_that("KrzCor returns correct results",
          {
              cov.matrix.1 <- RandCorr(10)
              cov.matrix.2 <- RandCorr(10)
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
            mat.list <- lapply(as.list(1:10), function(x) RandCorr(10))
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            results <- KrzCor(mat.list)
            results.2 <- KrzCor(mat.list, repeat.vector = rep.vec)
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
            lower.2  <- t(results.2)[upper.tri(results.2)]
            expect_that(upper.2, equals(upper))
            expect_that(diag(results.2), is_equivalent_to(rep.vec))
            expect_that(sum(lower.2 < upper.2), equals(0))
          }
)

test_that("KrzCor returns correct results on lists + matrices",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandCorr(10))
            y.matrix <- RandCorr(10)
            results <- KrzCor(mat.list, y.matrix)
            expect_that(results, is_a("numeric"))
            expect_that(sum(is.na(results)), equals(0))
            expect_that(length(results), equals(length(mat.list)))
            names(mat.list) <- 1:length(mat.list)
            named.results <- KrzCor(mat.list, y.matrix)
            expect_that(named.results, is_a("data.frame"))
            expect_that(sum(is.na(named.results)), equals(0))
            expect_that(dim(named.results), equals(c(length(mat.list), 2)))
            expect_that(named.results[,".id"], equals(names(mat.list)))
          }
)

