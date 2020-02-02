test_that("CalcEigenVar throws error",
          {
            expect_that(CalcEigenVar(array(1:100, c(10, 10))), 
                        throws_error("covariance matrix must be symmetric."))
          })
test_that("CalcEigenVar returns sensible results",
          {
            data(iris)
            cov.matrix = cov(iris[,1:4])
            eVals <- eigen(cov.matrix, only.values = TRUE)$values
            ident = diag(rep(2, 10))
            expect_that(CalcEigenVar(ident), equals(0))
            expect_that(CalcEigenVar(cov.matrix,sd = FALSE, rel=FALSE), 
                        equals(sum((eVals-mean(eVals))^2)/length(eVals)))
            expect_that(CalcEigenVar(cov.matrix,sd = TRUE, rel=FALSE), 
                        equals(sqrt(sum((eVals-mean(eVals))^2)/length(eVals))))
          }
)