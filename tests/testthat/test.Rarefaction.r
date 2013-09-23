test_that("Rarefaction returns sensible results",
          {
            ind.data <- matrix(rnorm(30*10), 30, 10)
            results <- Rarefaction(ind.data, num.reps = 5)
            expect_that(results, is_a("list"))
            expect_that(names(results), equals(as.character(1:29)))
            expect_that(length(results[[1]]), equals(5))
          }
)
