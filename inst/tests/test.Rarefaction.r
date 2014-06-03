test_that("RarefactionRandomSkewers returns sensible results",
          {
            ind.data <- matrix(rnorm(30*10), 30, 10)
            results <- Rarefaction(ind.data, "randomskewers", num.reps = 5)
            expect_that(results, is_a("list"))
            expect_that(names(results), equals(as.character(1:29)))
            expect_that(length(results[[1]]), equals(5))
          }
)
test_that("RarefactionMantelCor returns sensible results",
          {
            ind.data <- matrix(rnorm(30*10), 30, 10)
            results <- Rarefaction(ind.data, "mantel", num.reps = 5)
            expect_that(results, is_a("list"))
            expect_that(names(results), equals(as.character(1:29)))
            expect_that(length(results[[1]]), equals(5))
          }
)
test_that("RarefactionKrzCor returns sensible results",
          {
            ind.data <- matrix(rnorm(30*10), 30, 10)
            results <- Rarefaction(ind.data, "krz", num.reps = 5)
            expect_that(results, is_a("list"))
            expect_that(names(results), equals(as.character(1:29)))
            expect_that(length(results[[1]]), equals(5))
          }
)
test_that("RarefactionKrzCor returns sensible results",
          {
            ind.data <- matrix(rnorm(30*10), 30, 10)
            results <- Rarefaction(ind.data, "krz", num.reps = 5, correlation = TRUE)
            expect_that(results, is_a("list"))
            expect_that(names(results), equals(as.character(1:29)))
            expect_that(length(results[[1]]), equals(5))
          }
)
