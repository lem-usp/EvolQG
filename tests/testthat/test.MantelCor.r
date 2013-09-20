test_that("MantelCor returns correct results",
          {
            cor.matrix.1 <- cor(matrix(rnorm(30*10), 30, 10))
            cor.matrix.2 <- cor(matrix(rnorm(30*10), 30, 10))
            results <- MantelCor(cor.matrix.1, cor.matrix.2)
            expect_that(length(results), equals(2))
            expect_that(results[1] <=  1, is_true())
            expect_that(results[1] >= -1, is_true())
            expect_that(results[2] <=  1, is_true())
            expect_that(results[2] >=  0, is_true())
            expect_that(MantelCor(cor.matrix.1, cor.matrix.1), is_equivalent_to(c(1, 0.001)))
          }
)
