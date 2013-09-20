test_that("MultiKrzCor returns correct results",
          {
            cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.3 <- cov(matrix(rnorm(30*10), 30, 10))
            mat.list <- list(cov.matrix.1, cov.matrix.2, cov.matrix.3)
            results <- MultiKrzCor(mat.list)
            results.2 <- MultiKrzCor(mat.list, c(0.8, 0.9, 0.85))
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
            expect_that(diag(results.2), equals(c(0.8, 0.9, 0.85)))
            expect_that(sum(lower.2 > upper.2), equals(length(mat.list)))
          }
)

test_that("MultiKrzProjection returns correct results",
          {
            cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.3 <- cov(matrix(rnorm(30*10), 30, 10))
            mat.list <- list(cov.matrix.1, cov.matrix.2, cov.matrix.3)
            expect_that(dim(MultiKrzProjection(mat.list)), equals(c(6, 3)))
            expect_that(length(MultiKrzProjection(mat.list, full.results = T)), equals(3))
            expect_that(dim(MultiKrzProjection(mat.list, full.results = T)), equals(NULL))
          }
)

test_that("MultiMahalanobis returns correct results",
          {
            mean.1 <- colMeans(matrix(rnorm(30*10), 30, 10))
            mean.2 <- colMeans(matrix(rnorm(30*10), 30, 10))
            mean.3 <- colMeans(matrix(rnorm(30*10), 30, 10))
            mean.list <- list(mean.1, mean.2, mean.3)
            euclidian <- MultiMahalanobis(mean.list, diag(rep(1, 10)))
            half.euclidian <- MultiMahalanobis(mean.list, diag(rep(0.5, 10)))
            expect_that(euclidian*2, equals(half.euclidian))
            expect_that(euclidian, is_a("dist"))
            expect_that(length(euclidian), equals(3))
          }
)

test_that("MultiRSMantel returns corret results",
          {
            cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.3 <- cov(matrix(rnorm(30*10), 30, 10))
            mat.list <- list(cov.matrix.1, cov.matrix.2, cov.matrix.3)
            set.seed(42)
            results.list <- MultiRsMantel(mat.list)
            results <- results.list[[1]]
            probabilities <- results.list[[2]]
            expect_that(results.list, is_a("list"))
            set.seed(42)
            results.list.2 <- MultiRsMantel(mat.list, repeat.vector = c(0.8, 0.9, 0.85))
            results.2 <- results.list.2[[1]]
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
            expect_that(diag(results.2), equals(c(0.8, 0.9, 0.85)))
            expect_that(sum(lower.2 > upper.2), equals(length(mat.list)))
          }
)

test_that("MultiRSMantel returns corret results",
          {
            cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.3 <- cov(matrix(rnorm(30*10), 30, 10))
            mat.list <- list(cov.matrix.1, cov.matrix.2, cov.matrix.3)
            results.list <- MultiRsMantel(mat.list, MantelCor)
            results <- results.list[[1]]
            probabilities <- results.list[[2]]
            expect_that(results.list, is_a("list"))
            results.list.2 <- MultiRsMantel(mat.list, MantelCor, repeat.vector = c(0.8, 0.9, 0.85))
            results.2 <- results.list.2[[1]]
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
            expect_that(diag(results.2), equals(c(0.8, 0.9, 0.85)))
            expect_that(sum(abs(lower.2) > abs(upper.2)), equals(length(mat.list)))
          }
)
