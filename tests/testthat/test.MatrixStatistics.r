test_that("Autonomy returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              auto <- ((t (beta) %*% solve (cov.matrix) %*% beta)^(-1)) / (t (beta) %*% cov.matrix %*% beta)
              expect_that(Autonomy(cov.matrix, beta), equals(as.numeric(auto)))
              rmatrix <- RandomMatrix(40, 1, 1, 10)
              ev <- eigen(rmatrix)
              ev$values[35:40] <- 0
              singular.matrix <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
              expect_warning(Autonomy(singular.matrix), 
                             "matrix is singular, can't compute autonomy directly. Using nearPD, results could be wrong")
          }
)
test_that("ConditionalEvolvability returns correct restults",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              cond.evol <- (t (beta) %*% solve (cov.matrix) %*% beta)^(-1)
              expect_that(ConditionalEvolvability(cov.matrix, beta), equals(as.numeric(cond.evol)))
              rmatrix <- RandomMatrix(40, 1, 1, 10)
              ev <- eigen(rmatrix)
              ev$values[35:40] <- 0
              singular.matrix <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
              expect_warning(ConditionalEvolvability(singular.matrix), 
                             "matrix is singular, can't compute conditional evolvability directly. Using nearPD, results could be wrong")
          }
)
test_that("Constraints returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta = Normalize(rnorm(4))
              const = abs (t (Normalize (eigen (cov.matrix)$vectors[,1])) %*% Normalize (cov.matrix %*% beta))
              expect_that(Constraints(cov.matrix, beta), equals(as.numeric(const)))
          }
)
test_that("Evolvability returns correct results",
          {
              data(iris)
              cov.matrix <- cov(iris[,1:4])
              beta <- Normalize(rnorm(4))
              evol <- t (beta) %*% cov.matrix %*% beta
              expect_that(Evolvability(cov.matrix, beta), equals(as.numeric(evol)))
          }
)
test_that("Flexibility returns correct result",
         {
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           flex <- t (beta) %*% cov.matrix %*% beta / Norm (cov.matrix %*% beta)
           expect_that(Flexibility(cov.matrix, beta), equals(as.numeric(flex)))
           expect_true(Flexibility(cov.matrix, beta) <=  1)
           expect_true(Flexibility(cov.matrix, beta) >= -1)
         }
)

test_that("Pc1Percent returns correct results",
          {
            cov.matrix = cov(matrix(rnorm(30*10), 30, 10))
            pc1 <- eigen(cov.matrix)$values[1]/sum(diag(cov.matrix))
            expect_that(Pc1Percent(cov.matrix), equals(pc1))
            expect_true(Pc1Percent(cov.matrix) <=  1)
            expect_true(Pc1Percent(cov.matrix) > 0)
          }
)
test_that("Respondability returns correct result",
         {
           data(iris)
           beta <- Normalize(rnorm(4))
           cov.matrix <- cov(iris[,1:4])
           repond <- Norm(cov.matrix %*% beta)
           expect_that(Respondability(cov.matrix, beta), equals(repond))
         }
)

test_that("MeanMatrixStatistics excludes legacy integration metrics by default",
          {
            stats <- MeanMatrixStatistics(cov(iris[,1:4]))
            expect_false("MeanSquaredCorrelation" %in% names(stats))
            expect_false("ICV" %in% names(stats))
            expect_true("pc1.percent" %in% names(stats))
            expect_true("EigenSd" %in% names(stats))
          }
)

test_that("MeanMatrixStatistics can include R2 and ICV on demand",
          {
            cov.matrix <- cov(iris[,1:4])
            stats <- suppressWarnings(MeanMatrixStatistics(cov.matrix, include.R2 = TRUE, include.ICV = TRUE))
            expect_equal(as.numeric(stats["MeanSquaredCorrelation"]), mean(cov2cor(cov.matrix)[lower.tri(cov2cor(cov.matrix))]^2))
            expect_equal(as.numeric(stats["ICV"]), sd(eigen(cov.matrix, only.values = TRUE)$values) / mean(eigen(cov.matrix, only.values = TRUE)$values))
          }
)
