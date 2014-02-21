test_that("Ancestral Means gives correct results", 
         {
           TREE <- "iris((setosa:1,versicolor:1):1,virginica:2);"
           tree.iris <- read.tree(text = TREE)
           iris.cov.list <- dlply(iris, "Species", function(x) cov(x[,1:4]))
           cov.matrices <- AncestralStates(tree.iris, iris.cov.list)
           w_matrix <- CalculateMatrix(lm(as.matrix(iris[,1:4])~iris[,5]))
           expect_that(cov.matrices$'4', equals(w_matrix))
           versicolor_setosa <- Reduce("+", cov.matrices[1:2])/2
           expect_that(cov.matrices$'5', equals(versicolor_setosa))
           expect_that(AncestralStates(tree.iris, iris.cov.list[1:2]), 
                       throws_error("All tip labels must be in stat list."))
})