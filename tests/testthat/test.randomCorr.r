test_that("RandomCorr returns sensible results",
          {
            mat.list <- lapply(as.list(1:100), function(x) RandCorr(10))
            mean.corr <- mean(MantelCor(mat.list, iterations = 2)[[1]][upper.tri(array(0, c(100, 100)))])
            expect_that(dim(mat.list[[1]]), equals(c(10, 10)))
            expect_that(sum(diag(mat.list[[1]])), equals(10))
            expect_that(abs(mean.corr) < 0.01, is_true())
          }
)
