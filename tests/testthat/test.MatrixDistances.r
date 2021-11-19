test_that("Matrix Reimann Distance",
{
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(10)
  cor.matrix.1 <- RandomMatrix(10, LKJ = FALSE)
  cor.matrix.2 <- RandomMatrix(10, LKJ = FALSE)
  results <- RiemannDist(cor.matrix.1, cor.matrix.2)
  results_rec <- RiemannDist(cor.matrix.2, cor.matrix.1)
  results_ident <- RiemannDist(cor.matrix.1, cor.matrix.1)
  expect_that(results, is_a("numeric"))
  expect_true(results >=  0)
  expect_that(results, equals(results_rec, tolerance = 10e-4))
  expect_that(results_ident, equals(0))
}
)

test_that("Matrix Distribution Overlap Distance",
{
  cor.matrix.1 <- RandomMatrix(10, LKJ = FALSE)
  cor.matrix.2 <- RandomMatrix(10, LKJ = FALSE)
  results <- OverlapDist(cor.matrix.1, cor.matrix.2)
  results_rec <- OverlapDist(cor.matrix.2, cor.matrix.1)
  results_ident <- OverlapDist(cor.matrix.1, cor.matrix.1)
  expect_that(results, is_a("numeric"))
  expect_true(results >=  0)
  expect_that(results, equals(results_rec, tolerance = 10e-4))
  expect_that(results_ident, equals(0))
}
)

