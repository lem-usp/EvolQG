plot.group.evol <-
function (evol.list)
{
  upper.crit <- function (vec, val = .95)
  {
    n.obs <- length (vec)
    vec <- sort (vec)
    crit.pos <- round (val * n.obs)
    return (vec[crit.pos])
  }
  plot.single <- function (evol.list, which)
  {
    n.obj <- length (evol.list)
    extract.dist <- function (element, which) return (element$dist[,which]/element$max.val[which])
    to.plot <- sapply (evol.list, extract.dist, which)
    colnames (to.plot) <- names (evol.list)
    boxplot (to.plot, cex = 0.2, boxwex = 0.5, border = rgb (.6,.6,.6))
    extract.mod <- function (element, which) return (element$mod[,which]/element$max.val[which])
    extract.mod.names <- function (element) return (rownames (element$mod))
    mod.values <- lapply (evol.list, extract.mod, which)
    mod.names <-  lapply (evol.list, extract.mod.names)
    crits <- apply (to.plot, 2, upper.crit)
    for (i in 1:n.obj)
    {
      segments (x0 = i - 0.25, x1 = i + 0.25, y0 = crits [i],
                col = 'red', lwd = 2)
      segments (x0 = i - 0.25, x1 = i + 0.25, y0 = mod.values [[i]],
                col = 'black', lty = i)
    }
  }
}
