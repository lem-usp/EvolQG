test_that("AlphaRep returns correct values",
          {
              data(iris)
              cor.matrix = cor(iris[,1:4])
              tam = length(iris[,1])
              vec <- cor.matrix[lower.tri(cor.matrix)]
              var.erro <- (1 - mean(vec)^2)/(tam-2)
              var.vec <- var(vec)
              iris.rep = (var.vec - var.erro)/var.vec
              expect_that(AlphaRep(cor.matrix, tam), equals(iris.rep))
          }
)

