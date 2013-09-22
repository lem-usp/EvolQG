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
