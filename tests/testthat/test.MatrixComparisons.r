test_that("RandomSkewers returns correct results on pairs of matrices",
          {
            expect = c(correlation = 0.954184140494324, probability = 0.00561, correlation_sd = 0.0553934050116506)
            data(dentus)
            cor.matrix.1 <- cov(dentus[dentus$species=="A", 1:4])
            cor.matrix.2 <- cov(dentus[dentus$species=="B", 1:4])
            results <- RandomSkewers(cor.matrix.1, cor.matrix.2, num.vectors = 100000)
            expect_equal(results, expect, tolerance = 1e-2)
            expect_that(results, is_a("numeric"))
            expect_that(length(results), equals(3))
            expect_true(results[1] <=  1)
            expect_true(results[1] >= -1)
            expect_true(results[2] <=  1)
            expect_true(results[2] >=  0)
            expect_true(results[3] >=  0)
          }
)

test_that("RandomSkewers returns correct results on lists",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(10))
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            results.list <- RandomSkewers(mat.list)
            expect_that(results.list, is_a("list"))
            results <- results.list[[1]]
            expect_that(sum(is.na(results)), equals(0))
            probabilities <- results.list[[2]]
            expect_that(sum(is.na(probabilities)), equals(0))
            expect_that(results.list, is_a("list"))
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            results.list.2 <- RandomSkewers(mat.list, repeat.vector = rep.vec)
            results.2 <- results.list.2[[1]]
            expect_that(sum(is.na(results.2)), equals(0))
            probabilities.2 <- results.list.2[[2]]
            expect_that(dim(results), equals(c(length(mat.list),length(mat.list))))
            lower  <- results[lower.tri(results)]
            lower.bool = sapply(lower, function(x) isTRUE(x > -1 & x < 1))
            expect_that(sum(lower.bool), equals(length(lower.bool)))
            upper  <- results[upper.tri(results, diag = T)]
            upper.bool = sapply(upper, function(x) isTRUE(x == 0))
            expect_that(sum(upper.bool), equals(length(upper.bool)))
            lower.2  <- results.2[lower.tri(results.2)]
            upper.2  <- t(results.2)[lower.tri(results.2)]
            expect_that(lower.2, equals(lower))
            expect_that(diag(results.2), is_equivalent_to(rep.vec))
            expect_that(sum(upper.2 < lower.2), equals(0))
          }
)

test_that("RandomSkewers returns correct results on lists + matrices",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(10))
            y.matrix <- RandomMatrix(10)
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            results <- RandomSkewers(mat.list, y.matrix)
            expect_that(results, is_a("data.frame"))
            expect_that(sum(is.na(results)), equals(0))
            expect_that(dim(results), equals(c(length(mat.list), 3)))
            names(mat.list) <- 1:length(mat.list)
            named.results <- RandomSkewers(mat.list, y.matrix)
            expect_that(named.results, is_a("data.frame"))
            expect_that(sum(is.na(named.results)), equals(0))
            expect_that(dim(named.results), equals(c(length(mat.list), 4)))
            expect_that(named.results[,".id"], equals(names(mat.list)))
          }
)
test_that("MantelCor returns correct results on non-landmark data",
          {
            cor.matrix.1 <- RandomMatrix(15)
            cor.matrix.2 <- RandomMatrix(15)
            results <- MantelCor(cor.matrix.1, cor.matrix.2)
            expect_equivalent(results[1], cor(cor.matrix.1[lower.tri(diag(15))], 
                                              cor.matrix.2[lower.tri(diag(15))]))
            
          })

test_that("MantelCor returns correct results on landmark data",
          {
            cor.matrix.1 <- RandomMatrix(12)
            cor.matrix.2 <- RandomMatrix(12)
            results <- MantelCor(cor.matrix.1, cor.matrix.2, landmark.dim = 2)
            expect_equivalent(results[1], cor(cor.matrix.1[evolqg:::lower.tri.land(diag(12), 2)], 
                                              cor.matrix.2[evolqg:::lower.tri.land(diag(12), 2)]))
            results <- MantelCor(cor.matrix.1, cor.matrix.2, landmark.dim = 3)
            expect_equivalent(results[1], cor(cor.matrix.1[evolqg:::lower.tri.land(diag(12), 3)], 
                                              cor.matrix.2[evolqg:::lower.tri.land(diag(12), 3)]))
            expect_error(MantelCor(cor.matrix.1, cor.matrix.2, landmark.dim = 4))
          })

test_that("MantelCor returns correct results on list",
          {
            mat.list <- RandomMatrix(10, 10)
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            results.list <- MantelCor(mat.list)
            results <- results.list[[1]]
            expect_equivalent(results[2,1], MantelCor(mat.list[[1]], 
                                                      mat.list[[2]])[1])
            expect_equivalent(results[8,5], MantelCor(mat.list[[8]], 
                                                      mat.list[[5]])[1])
            expect_that(sum(is.na(results)), equals(0))
            probabilities <- results.list[[2]]
            expect_that(sum(is.na(probabilities)), equals(0))
            expect_that(results.list, is_a("list"))
            results.list.2 <- MantelCor(mat.list, repeat.vector = rep.vec)
            results.2 <- results.list.2[[1]]
            expect_that(sum(is.na(results.2)), equals(0))
            probabilities.2 <- results.list.2[[2]]
            expect_that(dim(results), equals(c(length(mat.list),length(mat.list))))
            lower <- results[lower.tri(results)]
            lower.bool = sapply(lower, function(x) isTRUE(x > -1 & x < 1))
            expect_that(sum(lower.bool), equals(length(lower.bool)))
            upper  <- results[upper.tri(results, diag = T)]
            upper.bool = sapply(upper, function(x) isTRUE(x == 0))
            expect_that(sum(upper.bool), equals(length(upper.bool)))
            lower.2 <- results.2[lower.tri(results.2)]
            upper.2 <- t(results.2)[lower.tri(results.2)]
            expect_that(lower.2, equals(lower))
            expect_that(diag(results.2), is_equivalent_to(rep.vec))
            expect_that(sum(abs(upper.2) < abs(lower.2)), equals(0))
          }
)

test_that("MantelCor returns correct results on list for landmark data",
          {
            mat.list <- RandomMatrix(12, 10)
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            results.list <- MantelCor(mat.list, landmark.dim = 2)
            results <- results.list[[1]]
            expect_equivalent(results[2,1], MantelCor(mat.list[[1]], 
                                                      mat.list[[2]], landmark.dim = 2)[1])
            expect_equivalent(results[8,5], MantelCor(mat.list[[8]], 
                                                      mat.list[[5]], landmark.dim = 2)[1])
            results.list <- MantelCor(mat.list, landmark.dim = 3)
            results <- results.list[[1]]
            expect_equivalent(results[2,1], MantelCor(mat.list[[1]], 
                                                      mat.list[[2]], landmark.dim = 3)[1])
            expect_equivalent(results[8,5], MantelCor(mat.list[[8]], 
                                                      mat.list[[5]], landmark.dim = 3)[1])
            expect_that(sum(is.na(results)), equals(0))
            probabilities <- results.list[[2]]
            expect_that(sum(is.na(probabilities)), equals(0))
            expect_that(results.list, is_a("list"))
            results.list.2 <- MantelCor(mat.list, repeat.vector = rep.vec, landmark.dim = 3)
            results.2 <- results.list.2[[1]]
            expect_that(sum(is.na(results.2)), equals(0))
            probabilities.2 <- results.list.2[[2]]
            expect_that(dim(results), equals(c(length(mat.list),length(mat.list))))
            lower <- results[lower.tri(results)]
            lower.bool = sapply(lower, function(x) isTRUE(x > -1 & x < 1))
            expect_that(sum(lower.bool), equals(length(lower.bool)))
            upper  <- results[upper.tri(results, diag = T)]
            upper.bool = sapply(upper, function(x) isTRUE(x == 0))
            expect_that(sum(upper.bool), equals(length(upper.bool)))
            lower.2 <- results.2[lower.tri(results.2)]
            upper.2 <- t(results.2)[lower.tri(results.2)]
            expect_that(lower.2, equals(lower))
            expect_that(diag(results.2), is_equivalent_to(rep.vec))
            expect_that(sum(abs(upper.2) < abs(lower.2)), equals(0))
          }
)

test_that("MantelCor returns correct results on lists + matrices",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(11))
            y.matrix <- RandomMatrix(11)
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            results <- MantelCor(mat.list, y.matrix)
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            expect_equivalent(results[1,], MantelCor(y.matrix, mat.list[[1]]))
            expect_that(results, is_a("data.frame"))
            expect_that(sum(is.na(results)), equals(0))
            expect_that(dim(results), equals(c(length(mat.list), 2)))
            names(mat.list) <- 1:length(mat.list)
            named.results <- MantelCor(mat.list, y.matrix)
            expect_that(named.results, is_a("data.frame"))
            expect_that(sum(is.na(named.results)), equals(0))
            expect_that(dim(named.results), equals(c(length(mat.list), 3)))
            expect_that(named.results[,".id"], equals(names(mat.list)))
          }
)

test_that("MantelCor returns correct results on lists + matrices for landmark data",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(12))
            y.matrix <- RandomMatrix(12)
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            results <- MantelCor(mat.list, y.matrix, landmark.dim = 2)
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            expect_equivalent(results[1,], MantelCor(y.matrix, mat.list[[1]], landmark.dim = 2))
            expect_that(results, is_a("data.frame"))
            expect_that(sum(is.na(results)), equals(0))
            expect_that(dim(results), equals(c(length(mat.list), 2)))
            names(mat.list) <- 1:length(mat.list)
            named.results <- MantelCor(mat.list, y.matrix)
            expect_that(named.results, is_a("data.frame"))
            expect_that(sum(is.na(named.results)), equals(0))
            expect_that(dim(named.results), equals(c(length(mat.list), 3)))
            expect_that(named.results[,".id"], equals(names(mat.list)))
          }
)

test_that("KrzCor returns correct results",
          {
              cov.matrix.1 <- RandomMatrix(10)
              cov.matrix.2 <- RandomMatrix(10)
              ret.dim = round(dim(cov.matrix.1)[1]/2 - 1)
              EigenVectors <- function (x) return (eigen(x)$vectors[,1:ret.dim])
              A <- EigenVectors (cov.matrix.1)
              B <- EigenVectors (cov.matrix.2)
              S <- t(A) %*% B %*% t(B) %*% A
              SL <- sum (eigen(S)$values) / ret.dim
              expect_that(KrzCor(cov.matrix.1, cov.matrix.2), is_equivalent_to(SL))
              expect_that(KrzCor(cov.matrix.1, cov.matrix.2, 10), is_equivalent_to(1))
              expect_that(KrzCor(x  <- RandomMatrix(11), x), is_equivalent_to(1))
              expect_true(KrzCor(cov.matrix.1, cov.matrix.2) <= 1)
              expect_true(KrzCor(cov.matrix.1, cov.matrix.2) > 0)
          }
)

test_that("KrzCor returns correct results",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(10))
            rep.vec <- runif(length(mat.list), 0.8, 0.9)
            results <- KrzCor(mat.list)
            expect_equivalent(results[2,1], KrzCor(mat.list[[1]], mat.list[[2]]))
            results.2 <- KrzCor(mat.list, repeat.vector = rep.vec)
            expect_that(sum(is.na(results)), equals(0))
            expect_that(sum(is.na(results.2)), equals(0))
            expect_that(dim(results), equals(c(length(mat.list),length(mat.list))))
            lower  <- results[lower.tri(results)]
            lower.bool = sapply(lower, function(x) isTRUE(x > 0 & x < 1))
            expect_that(sum(lower.bool), equals(length(lower.bool)))
            upper  <- results[upper.tri(results, diag = T)]
            upper.bool = sapply(upper, function(x) isTRUE(x == 0))
            expect_that(sum(upper.bool), equals(length(upper.bool)))
            lower.2  <- results.2[lower.tri(results.2)]
            upper.2  <- t(results.2)[lower.tri(results.2)]
            expect_that(lower.2, equals(lower))
            expect_that(diag(results.2), is_equivalent_to(rep.vec))
            expect_that(sum(upper.2 < lower.2), equals(0))
          }
)

test_that("KrzCor returns correct results on lists + matrices",
          {
            mat.list <- lapply(as.list(1:10), function(x) RandomMatrix(13))
            y.matrix <- RandomMatrix(13)
            results <- KrzCor(mat.list, y.matrix)
            expect_equivalent(results[1], KrzCor(mat.list[[1]], y.matrix))
            expect_that(results, is_a("numeric"))
            expect_that(sum(is.na(results)), equals(0))
            expect_that(length(results), equals(length(mat.list)))
            names(mat.list) <- 1:length(mat.list)
            named.results <- KrzCor(mat.list, y.matrix)
            expect_that(named.results, is_a("data.frame"))
            expect_that(sum(is.na(named.results)), equals(0))
            expect_that(dim(named.results), equals(c(length(mat.list), 2)))
            expect_that(named.results[,".id"], equals(names(mat.list)))
          }
)

test_that("KrzProjection returns correct results on matrices",
          {
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(17)
            cov.matrix.1 <- RandomMatrix(10, LKJ = FALSE)
            cov.matrix.2 <- RandomMatrix(10, LKJ = FALSE)
            expect_that(KrzProjection(cov.matrix.1, cov.matrix.2), is_a("list"))
            expect_that(length(KrzProjection(cov.matrix.1, cov.matrix.2)), equals(2))
            expect_that(length(KrzProjection(cov.matrix.1,
                                             cov.matrix.2,
                                             ret.dim.1 = 3)[[2]]),
                        equals(3))
            expect_that(length(KrzProjection(cov.matrix.1,
                                             cov.matrix.2,
                                             ret.dim.1 = 5)[[2]]),
                        equals(5))
            expect_that(length(KrzProjection(cov.matrix.1,
                                             cov.matrix.2,
                                             ret.dim.2 = 10)[[2]]),
                        equals(4))
            expect_that(KrzProjection(cov.matrix.1,
                                      cov.matrix.2,
                                      ret.dim.1 = 10, ret.dim = 10)[[1]],
                        equals(1))
            per.PC <- KrzProjection(cov.matrix.1, cov.matrix.2, ret.dim.1 = 10, ret.dim = 10)[[2]]
            ones = sapply(per.PC, function(x) isTRUE(all.equal(x, 1)))
            expect_that(sum(ones), equals(10))
          }
)

test_that("KrzProjection returns correct results on lists",
          {
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(170)
            mat.list <- RandomMatrix(10, 10, LKJ = FALSE)
            expect_that(dim(KrzProjection(mat.list)), equals(c(length(mat.list), length(mat.list))))
            expect_that(KrzProjection(mat.list)[1,2], equals(KrzProjection(mat.list[[1]], mat.list[[2]])[[1]]))
            expect_that(KrzProjection(mat.list)[2,1], equals(KrzProjection(mat.list[[2]], mat.list[[1]])[[1]]))
            expect_that(length(KrzProjection(mat.list, full.results = T)), equals(length(mat.list)))
            expect_that(KrzProjection(mat.list, full.results = T), is_a("list"))
          }
)

test_that("KrzProjection returns correct results on lists and matrices",
          {
            mat.list <- RandomMatrix(10, 10)
            results <- KrzProjection(mat.list, mat.list[[1]])
            test.results <- KrzProjection(mat.list)[,1]
            expect_that(results, is_a('data.frame'))
            expect_that(results[[1]], is_equivalent_to(test.results))
          }
)


