test_that("Matrix Reimann Distance",
{
  cor.matrix.1 <- RandomMatrix(10)
  cor.matrix.2 <- RandomMatrix(10)
  results <- MatrixRiemannDist(cor.matrix.1, cor.matrix.2)
  results_rec <- MatrixRiemannDist(cor.matrix.2, cor.matrix.1)
  results_ident <- MatrixRiemannDist(cor.matrix.1, cor.matrix.1)
  expect_that(results, is_a("numeric"))
  expect_that(results >=  0, is_true())
  expect_that(results, equals(results_rec))
  expect_that(results_ident, equals(0))
}
)

test_that("Matrix Distribution Overlap Distance",
{
  cor.matrix.1 <- RandomMatrix(10)
  cor.matrix.2 <- RandomMatrix(10)
  results <- MatrixOverlapDist(cor.matrix.1, cor.matrix.2)
  results_rec <- MatrixOverlapDist(cor.matrix.2, cor.matrix.1)
  results_ident <- MatrixOverlapDist(cor.matrix.1, cor.matrix.1)
  expect_that(results, is_a("numeric"))
  expect_that(results >=  0, is_true())
  expect_that(results, equals(results_rec))
  expect_that(results_ident, equals(0))
}
)