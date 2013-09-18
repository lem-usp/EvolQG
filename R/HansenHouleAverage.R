HansenHouleAverage <-
function (mat, num.vectors = 10000)
{
  with (load.of.functions,
        {
          n.char <- dim (mat) [1]
          beta.mat <- array (rnorm (n.char * num.vectors), c(n.char, num.vectors))
          beta.mat <- apply (beta.mat, 2, Normalize)
          iso.vec <- Normalize (rep(1, times = n.char))
          null.dist <- abs (t (iso.vec) %*% beta.mat)
          null.dist <- sort (null.dist)
          crit.value <- null.dist [round (0.95 * num.vectors)]
          cat ('critical value: ', crit.value, '\n')
          parm.dist <- array (0, c(num.vectors, 8))
          HansenHouleWrap <- function (hh.func) return (apply (beta.mat, 2, hh.func, cov.matrix = mat))
          parm.dist [,1:6] <- sapply (hansen.houle, HansenHouleWrap)
          parm.dist[,7] <- as.numeric (parm.dist[,5] > crit.value)
          parm.dist[,8] <- as.numeric (parm.dist[,6] > crit.value)
          parm.dist <- cbind (parm.dist, null.dist)
          colnames (parm.dist) <- c('resp',
                                    'evol',
                                    'cond.evol',
                                    'auto',
                                    'flex',
                                    'const',
                                    'flex.n',
                                    'const.n',
                                    'null.dist')
          parm.av <- colMeans (parm.dist)
          parm.av[7:8] <- parm.av[7:8] * num.vectors
          pc1 <- eigen (mat)$vectors[,1]
          HansenHouleWrapPc1 <- function (hh.func) return (hh.func (beta = pc1, cov.matrix = mat))
          maximum <- sapply (hansen.houle, HansenHouleWrapPc1)
          integration <- c (CalcR2 (mat), Pc1Percent (mat))
          names (integration) <- c ('MeanSquaredCorrelation', 'pc1%')
          parm.av <- c (integration, parm.av)
          return (list ('dist' = parm.dist, 'mean' = parm.av, 'max.val' = maximum))
        })
}
