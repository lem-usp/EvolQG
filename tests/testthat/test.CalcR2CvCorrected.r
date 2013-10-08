test_that("CalcR2CvCorrected returns sensible results on residuals",
          {
              results <- CalcR2CvCorrected(iris[,1:4], 100)
              expect_that(results, is_a('matrix'))
              expect_that(colnames(results), equals(c("r2", "eVals_cv", "mean_cv")))
              expect_that(sum(sapply(results[,1], function(x) !(x > 0 & x < 1))), equals(0))
              expect_that(dim(results), equals(c(100, 3)))
          })

test_that("CalcR2CvCorrected returns sensible results on models",
          {
              iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
              results <- CalcR2CvCorrected(iris.lm, 100)
              expect_that(results, is_a('matrix'))
              expect_that(colnames(results), equals(c("r2", "eVals_cv", "mean_cv")))
              expect_that(sum(sapply(results[,1], function(x) !(x > 0 & x < 1))), equals(0))
              expect_that(dim(results), equals(c(100, 3)))
          })
