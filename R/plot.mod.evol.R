plot.mod.evol <-
function (evo.out, new.dev = TRUE)
{
  require (MASS)
  with (evo.out,
        {
          n.hip <- nrow (mod)
          n.sk <- nrow (dist)
          if (new.dev)
            par (mfrow = c(2,3))
          out <- array (0, dim (mod))
          for (j in 1:6)
          {
            dist[,j] <- dist [,j] / max.val [j]
            mod [,j] <- mod [,j] / max.val [j]
            x.lim <- c (ifelse (min (mod[,j]) < min (dist[,j]), min (mod[,j]), min (dist[,j])),
                       ifelse (max (mod[,j]) > max (dist[,j]), max (mod[,j]), max (dist[,j])))
            truehist (dist[,j], prob = TRUE, border = 'grey', col = 'white', xlim = x.lim,
                      xlab = colnames (mod)[j], main = '')
            lines (density (dist[,j]))
            abline (v = sort (dist[,j]) [round (0.975 * n.sk)], lty = 2, col = 'red')
            abline (v = sort (dist[,j]) [round (0.025 * n.sk)], lty = 2, col = 'red')
            for (i in 1:n.hip)
            {
              abline (v = mod[i,j], lty = 3, col = rgb (.2,.2,.2))
              mtext (side = 3, at = mod[i,j], text = rownames (mod)[i], las = 3, cex = 0.7)
              out [i,j] <- sum (mod [i,j] < dist[,j])/n.sk
            }
          }
          dimnames (out) <- dimnames (mod)
          return (out)
        })
}
