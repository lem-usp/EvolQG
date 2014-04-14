test_that("RarefactionRandomSkewers returns sensible results",
          {
            ind.data <- matrix(rnorm(30*10), 30, 10)
            results <- RarefactionRandomSkewers(ind.data, num.reps = 5)
            expect_that(results, is_a("list"))
            expect_that(names(results), equals(as.character(1:29)))
            expect_that(length(results[[1]]), equals(5))
          }
)
test_that("RarefactionMantelCor returns sensible results",
          {
            ind.data <- matrix(rnorm(30*10), 30, 10)
            results <- RarefactionMantelCor(ind.data, num.reps = 5)
            expect_that(results, is_a("list"))
            expect_that(names(results), equals(as.character(1:29)))
            expect_that(length(results[[1]]), equals(5))
          }
)
test_that("RarefactionKrzCor returns sensible results",
          {
            ind.data <- matrix(rnorm(30*10), 30, 10)
            results <- RarefactionKrzCor(ind.data, num.reps = 5)
            expect_that(results, is_a("list"))
            expect_that(names(results), equals(as.character(1:29)))
            expect_that(length(results[[1]]), equals(5))
          }
)
test_that("RarefactionKrzCor returns sensible results",
          {
            ind.data <- matrix(rnorm(30*10), 30, 10)
            results <- RarefactionKrzCor(ind.data, T, num.reps = 5)
            expect_that(results, is_a("list"))
            expect_that(names(results), equals(as.character(1:29)))
            expect_that(length(results[[1]]), equals(5))
          }
)
