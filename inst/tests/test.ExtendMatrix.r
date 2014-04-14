test_that("ExtendMatrix returns correct results",
          {
              ind.data <- array(rnorm(100*11), c(100, 11))
              cov.matrix = cov(ind.data)
              cut.off = sample(2:11, 1)
              ext.matrix = ExtendMatrix(cov.matrix, cut.off)
              last.eval = eigen(cov.matrix)$values[cut.off]
              eVals = eigen(ext.matrix)$values
              extended = sapply(eVals, function(x) isTRUE(all.equal(x, last.eval)))
              expect_that(sum(extended), equals(11-cut.off+1))
              expect_that(sum(eVals>0), equals(length(eVals)))
              expect_that(ExtendMatrix(cov.matrix[1:9, 1:9]), throws_error("matrix is too small"))
          }
)
