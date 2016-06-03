test_that("KrzSubspace returns correct results",{
  kr.subspace_aguirre_bench <- function(Gs, vec){
    if (dim(Gs)[[1]] != dim(Gs)[[2]]){
      stop("G array must be of order n x n x m x MCMCsamp")
    }
    if (is.na(dim(Gs)[4])) {
      stop("There are no MCMCsamples")
    }
    n <- dim(Gs)[[1]]
    m <- dim(Gs)[[3]]
    MCMCsamp <- dim(Gs)[[4]] 
    if(length(vec) != m){stop("vec must have length = m")}
    h <- function (g, v){
      AA <- array(, c(n, n, m))  
      for (k in 1:m){
        g.vec <- eigen(g[,,k])$vectors[,1:(v[k])] 
        AA[,,k] <- g.vec %*% t(g.vec)
      }
      H <- apply(AA, 1:2, sum)
      list(H = H, AA = AA)
    }
    #internal function to calculate AA and H
    MCMC.H <- array(, c(n, n, MCMCsamp))
    dimnames(MCMC.H) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[4]])      
    MCMC.AA <- array(, c(n, n, m, MCMCsamp))
    dimnames(MCMC.AA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]], dimnames(Gs)[[4]])
    for (i in 1:MCMCsamp){
      kr <- h(Gs[,,,i], v = vec)
      MCMC.H[,,i] <- kr$H
      MCMC.AA[,,,i] <- kr$AA
    }	
    #calculate AA and H for the ith MCMC sample of the G array		
    avH <- apply(MCMC.H, 1:2, mean)
    rownames(avH) <- dimnames(Gs)[[1]]
    colnames(avH) <- dimnames(Gs)[[1]]
    #calculate the posterior mean H
    avAA <- apply(MCMC.AA, 1:3, mean)
    dimnames(avAA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]])
    #calculate the posterior mean AA
    avH.vec <- eigen(avH)$vectors
    #eigenanalysis of posterior mean H	
    proj<- function(a, b) t(b) %*% a %*% b
    #internal function to do projection
    avH.theta <- matrix(, n, m)
    for (i in 1:n){
      for (i in 1:n){
        avH.theta[i,] <- acos(sqrt(apply(avAA, 3, proj, b = avH.vec[,i]))) * (180/pi)
      }
    }
    #angles between the eigenvectors posterior mean H and the posterior mean subspaces of each population
    MCMC.H.val <- matrix(, MCMCsamp, n)
    colnames(MCMC.H.val) <- paste("h", 1:n, sep="")
    for (i in 1:n){
      MCMC.H.val[,i] <- apply(MCMC.H, 3, proj, b = avH.vec[,i])
    }
    #posterior distribution of the genetic variance for the eigenvectors of posterior mean H 
    MCMC.H.theta <- array(, c(n, m, MCMCsamp))
    rownames(MCMC.H.theta) <- paste("h", 1:n, sep="")
    colnames(MCMC.H.theta) <- dimnames(Gs)[[3]]
    for(i in 1:n){
      for(j in 1:MCMCsamp){
        MCMC.H.theta[i,,j] <- acos(sqrt(apply(MCMC.AA[,,,j], 3, proj, b = avH.vec[,i]))) * (180/pi)
      }
    }
    #posterior distribution of the angles between the eigenvectors of posterior mean H and the MCMC samples of the subspaces of each population
    list(avAA = avAA, avH = avH, MCMC.AA = MCMC.AA, MCMC.H = MCMC.H, MCMC.H.val = MCMC.H.val, MCMC.H.theta = MCMC.H.theta)
  }
  cov.matrices = aperm(aaply(1:10, 1, function(x) laply(RandomMatrix(6, 40, variance = runif(6,1, 10)), 
                                                        identity)), c(3, 4, 1, 2))
  Hs = llply(alply(cov.matrices, 4, function(x) alply(x, 3)), function(x) KrzSubspace(x, 3)$H)
  avgH = Reduce("+", Hs)/length(Hs)
  avgH.vec <- eigen(avgH)$vectors
  MCMC.H.val = laply(Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))
  x = kr.subspace_aguirre_bench(cov.matrices, vec = rep(3, 10))
  expect_that(HPDinterval(as.mcmc(x$MCMC.H.val)), is_equivalent_to(HPDinterval(as.mcmc(MCMC.H.val))))
})