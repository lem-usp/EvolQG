test_that("RSProjection returns correct results",
          {
            rs_benchmark_aguirre <- function(Gs, p = 0.8, vec = 1000){
              if (dim(Gs)[[1]] != dim(Gs)[[2]]){
                stop("G array must be of order n x n x m x MCMCsamp")
              }
              if (is.na(dim(Gs)[4])) {
                stop("There are no MCMCsamples")
              }
              n <- dim(Gs)[[1]]
              m <- dim(Gs)[[3]]
              MCMCsamp <- dim(Gs)[[4]]
              rand.vec <- aaply(matrix(rnorm(n*vec), n, vec), 2, Normalize)
              #generate unit length random vectors  
              proj<- function(G,b) t(b) %*% G %*% (b)
              #internal function to do projection
              G.proj <- array(,c(MCMCsamp, m, vec))
              colnames(G.proj) <- dimnames(Gs)[[3]]
              for (i in 1:vec){
                G.proj[,,i]<- t(apply(Gs, 3:4, proj, b = rand.vec[i,]))
              }
              #project each random vector through each MCMC sample of each G
              prs <- cbind(rep(1:m, each = m), 1:m) 
              prs.comp <- prs[prs[,1] < prs[,2], , drop = FALSE] 
              #setting up an index for HPD comparisons
              proj.score <-matrix(,vec,((m^2 - m)/2))
              for (k in 1:vec){
                HPD.int <- HPDinterval(as.mcmc(G.proj[,,k]), prob = p)
                proj.score[k,] <- ifelse(HPD.int[prs.comp[,1],1] > HPD.int[prs.comp[,2],2] | HPD.int[prs.comp[,2],1] > HPD.int[prs.comp[,1],2],1,0) 
              }
              #for a given random vector, examine if the HPD intervals of any pair of G matrices overlap
              vec.score <-cbind(rand.vec, proj.score)
              colnames(vec.score) <- c(1:n, paste(dimnames(Gs)[[3]][prs.comp[, 1]], ".vs.", dimnames(Gs)[[3]][prs.comp[, 2]], sep = ""))
              #collate the random vectors and the outcome of their projection on the G matrices
              sig.vec <- subset(vec.score, rowSums(vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0) 
              #collate just the random vectors that resulted in significant differences in variance
              if(dim(sig.vec)[1] <= 1) {warning("There were <= 1 significant vectors, try a larger vec or lower p"); eig.R <- "Na"}
              else{
                eig.R <- eigen(cov(sig.vec[,1:n]))
                rownames(eig.R$vectors) <- dimnames(Gs)[[1]]
                colnames(eig.R$vectors) <- c(paste("e", 1:n, sep = ""))
              }  
              #eigen analysis of the R matrix
              list(G.proj = G.proj, vec.score = vec.score, eig.R = eig.R)
            }
            
            cov.matrices = aperm(aaply(1:10, 1, function(x) 
              laply(RandomMatrix(6, 40, 
                                 variance = runif(6, 1, 10)), 
                    identity)), 
              c(3, 4, 1, 2))
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            rs_proj = RSProjection(cov.matrices, p = 0.8)  
            suppressWarnings(RNGversion("3.5.0"))
            set.seed(42)
            rs_proj_bench = rs_benchmark_aguirre(cov.matrices, p = 0.8) 
            expect_that(rs_proj, is_equivalent_to(rs_proj_bench))
            }
)
  
            
            