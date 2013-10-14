test_that("TestModularity returns correct results",
          {
            cor.matrix <- RandomMatrix(10)
            rand.hipots <- matrix(sample(c(1, 0), 30, replace=T), 10, 3)
            hip.array <- Morphometrics::CreateHipotMatrix(rand.hipots)
            expect_that(sum(laply(hip.array, isSymmetric)), equals(4))
            expect_that(hip.array, is_a("list"))
            expect_that(dim(hip.array[[1]]), equals(c(10, 10)))
            mod.test <- TestModularity(cor.matrix, rand.hipots)
            expect_that(dim(mod.test), equals(c(4, 6)))
            expect_that(mod.test[,".id"], equals(c("1", "2", "3", "Full Integration")))
            expect_that(colnames(mod.test), equals(c(".id", "Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio")))
            expect_that(sum(mod.test[1:4, 2:5] >= -1), equals(16))
            expect_that(sum(mod.test[1:4, 2:5] <= 1), equals(16))
            expect_that(mod.test[,4]/mod.test[,5], equals(mod.test[,6]))
          }
)
