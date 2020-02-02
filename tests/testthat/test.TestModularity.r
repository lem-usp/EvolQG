test_that("TestModularity returns correct results",
          {
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(43)
            cor.matrix <- RandomMatrix(10)
            rand.hipots <- matrix(sample(c(1, 0), 30, replace=T), 10, 3)
            hip.array <- CreateHypotMatrix(rand.hipots)
            expect_that(sum(laply(hip.array, isSymmetric)), equals(4))
            expect_that(hip.array, is_a("list"))
            expect_that(dim(hip.array[[1]]), equals(c(10, 10)))
            mod.test <- TestModularity(cor.matrix, rand.hipots)
            expect_that(dim(mod.test), equals(c(4, 6)))
            expect_that(mod.test[,"hypothesis"], equals(c("1", "2", "3", "Full Integration")))
            expect_that(colnames(mod.test), equals(c("hypothesis", "Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio")))
            expect_true(all(mod.test[1:4, 2:5] >= -1))
            expect_true(all(mod.test[1:4, 2:5] <= 1))
            expect_that(mod.test[,4]/mod.test[,5], equals(mod.test[,6]))
          }
)

test_that("TestModularity returns correct results for Modularity Hypothesis Index",
          {
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(43)
            cov.matrix <- RemoveSize(RandomMatrix(11))
            rand.hipots <- matrix(sample(c(1, 0), 33, replace=T), 11, 3)
            hip.array <- CreateHypotMatrix(rand.hipots)
            mod.test <- TestModularity(cov.matrix, rand.hipots, MHI = TRUE)
            expect_that(dim(mod.test), equals(c(4, 6)))
            expect_that(mod.test[,"hypothesis"], equals(c("1", "2", "3", "Full Integration")))
            expect_that(colnames(mod.test), equals(c("hypothesis", "Rsquared", 
                                                     "Probability", "AVG+", "AVG-", "MHI")))
            expect_true(all(mod.test[1:4, 2:5] >= -1))
            expect_true(all(mod.test[1:4, 2:5] <= 1))
            expect_that((mod.test[,4] - mod.test[,5])/CalcEigenVar(cov.matrix,sd = TRUE, rel = FALSE), equals(mod.test[,6]))
          }
)

test_that("MantelModTest returns correct results for Modularity Hypothesis Index",
          {
            # Create a single modularity hypothesis:
            hypot = rep(c(1, 0), each = 6)
            cor.hypot = CreateHypotMatrix(hypot)
            
            # First with an unstructured matrix:
            expect = c(Rsquared = -0.0393920752173473, Probability = 0.556443556443556, 
                       `AVG+` = -0.0596580385711816, `AVG-` = -0.0271784073214974, MHI = -0.02718047)
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            un.cor = RandomMatrix(12, LKJ = FALSE)
            result = MantelModTest(cor.hypot, RemoveSize(un.cor), MHI = TRUE)
            expect_equal(result, expect)
            
            # Now with a modular matrix:
            expect = c(Rsquared = 1, Probability = 0.001998001998002, `AVG+` = 0.8, 
                       `AVG-` = 0.3, MHI = 0.325128044)
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            hypot.mask = matrix(as.logical(cor.hypot), 12, 12)
            mod.cor = matrix(NA, 12, 12)
            mod.cor[ hypot.mask] = 0.8 # within-modules
            mod.cor[!hypot.mask] = 0.3 # between-modules
            diag(mod.cor) = 1
            result = MantelModTest(cor.hypot, mod.cor, MHI = TRUE)
            expect_equal(result, expect)
            avg = CalcAVG(cor.hypot, mod.cor, MHI = TRUE)
            expect_equal(avg, result[3:5])
          }
)

test_that("MantelModTest returns correct results for non-landmark data",
          {
            # Create a single modularity hypothesis:
            hypot = rep(c(1, 0), each = 6)
            cor.hypot = CreateHypotMatrix(hypot)
             
            # First with an unstructured matrix:
            expect = structure(c(0.0810901974853795, 0.275724275724276, 0.158055575956429, 0.061087664745461, 2.587356655636), .Names = c("Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio")) 
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            un.cor = RandomMatrix(12, LKJ = FALSE)
            result = MantelModTest(cor.hypot, un.cor)
            expect_equal(result, expect)
            
            # Now with a modular matrix:
            expect = structure(c(0.99102012957595, 0.001998001998002, 0.87040471785857, 0.354135870234584, 2.45782704045767), .Names = c("Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio"))
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            hypot.mask = matrix(as.logical(cor.hypot), 12, 12)
            mod.cor = matrix(NA, 12, 12)
            mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
            mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
            diag(mod.cor) = 1
            
            result = MantelModTest(cor.hypot, mod.cor)
            expect_equal(result, expect)
            mod.cor<-(mod.cor+t(mod.cor))/2
            result = MantelModTest(cor.hypot, mod.cor, MHI = TRUE)
            expect_equal(CalcAVG(cor.hypot, mod.cor), result[3:5])
          })

test_that("MantelModTest returns correct results for landmark data",
          {
            # Create a single modularity hypothesis:
            hypot = rep(c(1, 0), each = 6)
            cor.hypot = CreateHypotMatrix(hypot)
            
            # First with an unstructured matrix:
            expect = structure(c(0.112390456689369, 0.259, 0.209791381269483, 0.0703378394237725, 2.98262475771438, -0.0679013877016532), .Names = c("Rsquared", "Probability", "AVG+", "AVG-", "AVG Ratio", "AVG within landmark"))
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            un.cor = RandomMatrix(12, LKJ = FALSE)
            result = MantelModTest(cor.hypot, un.cor, landmark.dim = 2)
            expect_equal(result, expect)
            
            # Now with a modular matrix:
            expect = structure(c(1, 0, 0.8, 0.3, 2.66666666666667, 0.1), .Names = c("Rsquared","Probability", "AVG+", "AVG-", "AVG Ratio", "AVG within landmark"))
            hypot.mask = matrix(as.logical(cor.hypot), 12, 12)
            mod.cor = matrix(NA, 12, 12)
            mod.cor[ hypot.mask] = 0.8 # within-modules
            mod.cor[!hypot.mask] = 0.3 # between-modules
            mod.cor[evolqg:::CreateWithinLandMat(6, 2)] = 0.1# within-land
            diag(mod.cor) = 1
            result = MantelModTest(cor.hypot, mod.cor, landmark.dim = 2)
            expect_equal(result, expect)
            expect_equal(CalcAVG(cor.hypot, mod.cor, landmark.dim = 2, MHI = FALSE), result[3:6])
            
            mod.cor = matrix(NA, 12, 12)
            mod.cor[ hypot.mask] = 0.8 # within-modules
            mod.cor[!hypot.mask] = 0.3 # between-modules
            mod.cor[evolqg:::CreateWithinLandMat(4, 3)] = 0.1 # within-land
            diag(mod.cor) = 1
            result = MantelModTest(cor.hypot, mod.cor, landmark.dim = 3)
            expect_equal(result, expect)
            expect_equal(CalcAVG(cor.hypot, mod.cor, landmark.dim = 3, MHI = FALSE), result[3:6])
            })

test_that("MantelModTest trows errors",
          {
            hypot = rep(c(1, 0), each = 6)
            cor.hypot = CreateHypotMatrix(hypot)
            un.cor = RandomMatrix(12)
            expect_error(MantelModTest(un.cor, cor.hypot), 
                         "modularity hypothesis matrix should be binary")
            expect_error(MantelModTest(cor.hypot, cor.hypot, landmark.dim = 4), 
                         "landmark.dim should be either 2 or 3 dimensions")
            expect_error(CalcAVG(un.cor, cor.hypot), 
                         "modularity hypothesis matrix should be binary")
            expect_error(CalcAVG(cor.hypot, cor.hypot, landmark.dim = 4), 
                         "landmark.dim should be either 2 or 3 dimensions")
          })