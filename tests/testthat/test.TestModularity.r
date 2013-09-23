test_that("TestModularity returns correct results",
          {
            cor.matrix <- cor(matrix(rnorm(30*10), 30, 10))
            rand.hipots <- matrix(sample(c(1, 0), 30, replace=T), 10, 3)
            hip.array <- CreateHipotMatrix(rand.hipots)
            expect_that(sum(apply(hip.array, 3, isSymmetric)), equals(4))
            expect_that(dim(hip.array), equals(c(10, 10, 4)))
            mod.test <- TestModularity(cor.matrix, rand.hipots)
            expect_that(dim(mod.test), equals(c(5, 4)))
            expect_that(colnames(mod.test), equals(c("1", "2", "3", "Full Integration")))
            expect_that(rownames(mod.test), equals(c("Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio")))
            expect_that(sum(mod.test[1:4, 1:4] >= -1), equals(16))
            expect_that(sum(mod.test[1:4, 1:4] <= 1), equals(16))
            expect_that(mod.test[3,]/mod.test[4,], equals(mod.test[5,]))
          }
)
