test_that("Autonomy returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              auto <- ((t (beta) %*% solve (cov.matrix) %*% beta)^(-1)) / (t (beta) %*% cov.matrix %*% beta)
              expect_that(Autonomy(beta, cov.matrix), equals(auto))
          }
)
test_that("ConditionalEvolvability returns correct restults",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              cond.evol <- (t (beta) %*% solve (cov.matrix) %*% beta)^(-1)
              expect_that(ConditionalEvolvability(beta, cov.matrix), equals(cond.evol))
          }
)
test_that("Constraints returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta = Normalize(rnorm(4))
              const = abs (t (Normalize (eigen (cov.matrix)$vectors[,1])) %*% Normalize (cov.matrix %*% beta))
              expect_that(Constraints(beta, cov.matrix), equals(const))
          }
)
test_that("Evolvability returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              evol <- t (beta) %*% cov.matrix %*% beta
              expect_that(Evolvability(beta, cov.matrix), equals(evol))
          }
)
test_that("Flexibility returns correct result",
         {
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           flex <- t (beta) %*% cov.matrix %*% beta / Norm (cov.matrix %*% beta)
           expect_that(Flexibility(beta, cov.matrix), equals(flex))
           expect_that(Flexibility(beta, cov.matrix) <=  1, is_true())
           expect_that(Flexibility(beta, cov.matrix) >= -1, is_true())
         }
)

test_that("Pc1Percent returns correct results",
          {
            cov.matrix = cov(matrix(rnorm(30*10), 30, 10))
            pc1 <- eigen(cov.matrix)$values[1]/sum(diag(cov.matrix))
            expect_that(Pc1Percent(cov.matrix), equals(pc1))
            expect_that(Pc1Percent(cov.matrix) <=  1, is_true())
            expect_that(Pc1Percent(cov.matrix) > 0, is_true())
          }
)
test_that("Respondability returns correct result",
         {
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           repond <- Norm(cov.matrix %*% beta)
           expect_that(Respondability(beta, cov.matrix), equals(repond))
         }
)

test_that("MeanMatrixStatistics returns correct results",
          {
            set.seed(42)
            iris.stats <- MeanMatrixStatistics(cov(iris[,1:4]))
            test.values <- read.table("iris.stats")
            expect_that(iris.stats, is_equivalent_to(test.values[,1]))
          }
)
