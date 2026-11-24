test_that("DriftTest returns resonable results",
{
  skip_on_cran()
  means <- array(rnorm(40*10), c(10, 40)) 
  means.list <- alply(means, 1)
  cov.matrix <- RandomMatrix(40, 1, 1, 10, LKJ = FALSE)
  test.array <- DriftTest(means, cov.matrix, FALSE)
  test.list <- DriftTest(means.list, cov.matrix, FALSE)
  expect_equal(test.array[-6], test.list[-6])
  expect_is(test.array, "list")
  expect_is(test.array[[1]], "lm")
  expect_is(test.array[[2]], "matrix")
  expect_is(test.array[[3]], "numeric")
  expect_is(test.array[[4]], "numeric")
  expect_is(test.array[[5]], "logical")
  expect_is(test.array[[6]], c("gg", "ggplot"))
  expect_equal(names(test.array), c("regression",
                                    "coefficient_CI_95",
                                    "log.between_group_variance",
                                    "log.W_eVals",
                                    "drift_rejected",
                                    "plot"))
  expect_equal(length(test.array$log.between_group_variance), dim(means)[2])
  expect_equal(length(test.array$log.W_eVals), dim(means)[2])
  expect_equal(test.array$log.W_eVals, log(eigen(cov.matrix)$values))
})

test_that("TreeDriftTest returns resonable results",
{
  library(ape)
  data(bird.orders)
  tree <- bird.orders
  mean.list <- llply(tree$tip.label, function(x) rnorm(5))
  names(mean.list) <- tree$tip.label
  cov.matrix.list <- RandomMatrix(5, length(tree$tip.label))
  names(cov.matrix.list) <- tree$tip.label
  w.cov <- PhyloW(tree, cov.matrix.list)$'24'
  test.list <- TreeDriftTest(tree, mean.list, cov.matrix.list)
  expect_is(test.list, 'list')
  expect_equal(length(test.list), 13)
  expect_equal(test.list[[length(test.list)]][-6], 
               DriftTest(mean.list, w.cov, FALSE)[-6])
})


# gets descendant node numbers
# written by Liam Revell 2012, 2013, 2014
getDescendants<-function(tree,node,curr=NULL){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  if(length(curr)==0&&node<=length(tree$tip.label)) curr<-node
  w<-which(daughters>length(tree$tip.label))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

test_that("TreeDriftTest returns resonable results",
          {
            library(ape)
            data(bird.orders)
            tree <- bird.orders
            mean.list <- llply(tree$tip.label, function(x) rnorm(5))
            names(mean.list) <- tree$tip.label
            cov.matrix.list <- RandomMatrix(5, length(tree$tip.label))
            names(cov.matrix.list) <- tree$tip.label
            cov.matrix.list = cov.matrix.list[sample(1:length(cov.matrix.list))]
            mean.list_2 = mean.list[sample(1:length(mean.list))]
            w.cov <- PhyloW(tree, cov.matrix.list)$'32'
            test.list <- TreeDriftTest(tree, mean.list_2, cov.matrix.list)
            expect_is(test.list, 'list')
            expect_equal(length(test.list), 13)
            expect_equal(test.list[["32"]][-6], 
                         DriftTest(mean.list[na.omit(tree$tip.label[getDescendants(tree, "32")])], 
                                   w.cov, FALSE)[-6])
          })

          
test_that("DriftTest sampling correction works",
{
  set.seed(42)
  p <- 6; m <- 21
  cov.matrix <- RandomMatrix(p, 1, 1, 10, LKJ = FALSE)
  means <- array(rnorm(m * p, sd = 3), c(m, p))
  n <- rep(30, m)
  uncorrected <- DriftTest(means, cov.matrix, FALSE)
  expect_equal(uncorrected[-6], DriftTest(means, cov.matrix, FALSE, sample.sizes = NULL)[-6])
  corrected <- DriftTest(means, cov.matrix, FALSE, sample.sizes = n)
  expect_equal(names(corrected), names(uncorrected))
  expected <- exp(uncorrected$log.between_group_variance) - eigen(cov.matrix)$values / 30
  expected[expected <= 0] <- NA
  expect_equal(corrected$log.between_group_variance, log(expected))
  expect_error(DriftTest(means, cov.matrix, FALSE, sample.sizes = n[-1]))
  # pure sampling error: most axes have no divergence beyond sampling and are dropped
  zero.means <- t(replicate(m, colMeans(MASS::mvrnorm(30, rep(0, p), cov.matrix))))
  expect_warning(try(DriftTest(zero.means, cov.matrix, FALSE, sample.sizes = n), silent = TRUE))
})
