test_that("MantelCor returns correct results",
          {
            cor.matrix.1 <- cor(matrix(rnorm(30*10), 30, 10))
            cor.matrix.2 <- cor(matrix(rnorm(30*10), 30, 10))
            results <- MantelCor(cor.matrix.1, cor.matrix.2)
            expect_that(length(results), equals(2))
            expect_that(results[1] <=  1, is_true())
            expect_that(results[1] >= -1, is_true())
            expect_that(results[2] <=  1, is_true())
            expect_that(results[2] >=  0, is_true())
          }
)
test_that("MantelCor returns corret results",
          {
            cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.3 <- cov(matrix(rnorm(30*10), 30, 10))
            mat.list <- list(cov.matrix.1, cov.matrix.2, cov.matrix.3)
            results.list <- MantelCor(mat.list)
            results <- results.list[[1]]
            expect_that(sum(is.na(results)), equals(0))
            probabilities <- results.list[[2]]
            expect_that(sum(is.na(probabilities)), equals(0))
            expect_that(results.list, is_a("list"))
            results.list.2 <- MantelCor(mat.list, repeat.vector = c(0.8, 0.9, 0.85))
            results.2 <- results.list.2[[1]]
            expect_that(sum(is.na(results.2)), equals(0))
            probabilities.2 <- results.list.2[[2]]
            expect_that(dim(results), equals(c(length(mat.list),length(mat.list))))
            upper  <- results[upper.tri(results)]
            upper.bool = sapply(upper, function(x) isTRUE(x > -1 & x < 1))
            expect_that(sum(upper.bool), equals(length(upper.bool)))
            lower  <- results[lower.tri(results, diag = T)]
            lower.bool = sapply(lower, function(x) isTRUE(x == 0))
            expect_that(sum(lower.bool), equals(length(lower.bool)))
            upper.2  <- results.2[upper.tri(results.2)]
            lower.2  <- results.2[lower.tri(results.2)]
            expect_that(upper.2, equals(upper))
            expect_that(diag(results.2), is_equivalent_to(c(0.8, 0.9, 0.85)))
            expect_that(sum(abs(lower.2) > abs(upper.2)), equals(length(mat.list)))
          }
)

