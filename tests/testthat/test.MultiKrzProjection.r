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
