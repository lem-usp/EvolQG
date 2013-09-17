test_that("BootstrapRep returns sensible results",
          {
              set.seed(42)
              ind.data <- iris[,1:4]
              nb <- 50
              n.ind <-  dim (ind.data) [1]
              original.cov.matrix <- var (ind.data)
              v.rep <- c()
              for (N in 1:nb){
                  sampled.data <- sample (1:n.ind, n.ind, TRUE)
                  sampled.data.cov.matrix <- var (ind.data[sampled.data,])
                  v.rep [N] <- RandomSkewers (original.cov.matrix, sampled.data.cov.matrix, 1000) [1]
              }
              out <- mean (v.rep)
              set.seed(42)
              expect_that(BootstrapRep(ind.data, 50), equals(out))
          }
)
