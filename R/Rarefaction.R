#' Rarefaction analysis via ressampling
#' Calculates the repeatability of a statistic of the data, such as
#' correlation or covariance matrix, via bootstrap ressampling with
#' varying sample sizes, from 2 to the size of the original data.
#'
#' Samples of various sizes, with replacement, are taken from the full population, a statistic calculated
#' and compared to the full population statistic. Prepackaged functions
#' for common comparison functions and statistics are suplied.
#'
#' A specialized ploting function displays the results in publication quality.
#' @aliases  Rarefaction RarefactionRandomSkewers RarefactionMantelCor RarefactionKrzCor PlotRarefaction
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param num.reps number of populations sampled per sample size
#' @param iterations Parameter for comparison function. Number of random skewers or number of permutations in mantel
#' @param ComparisonFunc Comparison function for calculated statistic
#' @param StatFunc Statistic to be calculated
#' @param num.cores Number of threads to use in computation. Requires doMC library.
#' @param correlation When using RarefactionKrzCor, statistic can be correlation or covariance. If TRUE, uses correlation.
#' @param ret.dim When using RarefactionKrzCor, number o retained dimensions may be specified
#' @param comparison.list output from rarefaction functions can be used in ploting
#' @param y.axis Y axis lable in plot
#' @details Bootstraping may be misleading with very small sample sizes. Use with caution.
#' @return returns the mean value of comparisons from samples to original statistic, for all sample sizes.
#' @author Diogo Melo, Guilherme Garcia
#' @seealso \code{\link{BootstrapRep}}
#' @rdname Rarefaction
#' @export
#' @import plyr
#' @importFrom ggplot2 ggplot aes layer scale_x_continuous scale_y_continuous theme_bw
#' @importFrom reshape2 melt
#' @examples
#' ind.data <- matrix(rnorm(30*10), 30, 10)
#' comparison.list <- Rarefaction(ind.data,
#'                                StatFunc = cov,
#'                                ComparisonFunc = function(x, y) RandomSkewers(x, y, 1000)[1],
#'                                num.reps=5,
#'                                num.cores = 1)
#'
#' results.RS <- RarefactionRandomSkewers(ind.data, num.reps = 5)
#' results.Mantel <- RarefactionMantelCor(ind.data, num.reps = 5)
#' results.KrzCov <- RarefactionKrzCor(ind.data, num.reps = 5)
#' results.KrzCor <- RarefactionKrzCor(ind.data, TRUE, num.reps = 5)
#'
#' #Easy access
#' library(reshape2)
#' melt(results.RS)
#'
#' #Plotting using ggplot2
#' a <- PlotRarefaction(results.RS, "Random Skewers")
#' b <- PlotRarefaction(results.Mantel, "Mantel")
#' c <- PlotRarefaction(results.KrzCov, "Krz - Covariance")
#' d <- PlotRarefaction(results.KrzCor, "Krz - Correlation")
#'
#' library(grid)
#' grid.newpage()
#' pushViewport(viewport(layout = grid.layout(2, 2)))
#' vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#' print(a, vp = vplayout(1, 1))
#' print(b, vp = vplayout(1, 2))
#' print(c, vp = vplayout(2, 1))
#' print(d, vp = vplayout(2, 2))
#'
#' @keywords rarefaction
#' @keywords bootstap
#' @keywords repeatability

Rarefaction <- function(ind.data,
                        StatFunc,
                        ComparisonFunc,
                        num.reps = 10,
                        num.cores = 1)
{
  if (num.cores > 1) {
    library(doMC)
    library(foreach)
    registerDoMC(num.cores)
    parallel = TRUE
  }
  else{
    parallel = FALSE
  }
  observed.stat = StatFunc(ind.data)
  num.ind = dim(ind.data)[1]
  MapStatFunc <- function(n){
    SampleFunction <- function(x){
      local.sample = sample(1:num.ind, n, replace=T)
      out <- StatFunc(ind.data[local.sample,])
      return(out)
    }
    return(alply(1:num.reps, 1, SampleFunction))
  }
  sample.stats = alply(2:num.ind, 1, MapStatFunc, .parallel = parallel)
  MapComparisonFunc <- function(stat.list){
    laply(stat.list, function(x) ComparisonFunc(x, observed.stat))
  }
  comparison.list = llply(sample.stats, MapComparisonFunc, .parallel = parallel)
  return(comparison.list)
}

#' @rdname Rarefaction
#' @export
PlotRarefaction <- function(comparison.list, y.axis = "Statistic"){
  plot.df = melt(comparison.list, value.name = 'avg.corr')
  plot.df = as.data.frame(lapply(plot.df, as.numeric))
  rarefaction.plot = ggplot(plot.df, aes(L1, avg.corr, group = L1)) +
  layer(geom = "boxplot") +
  scale_x_continuous("ReSample Size") +
  scale_y_continuous(y.axis) +
  theme_bw()
  return(rarefaction.plot)
}

#' @rdname Rarefaction
#' @export
RarefactionRandomSkewers <- function(ind.data,
                                     iterations = 1000,
                                     num.reps = 10,
                                     num.cores = 1){
  comparison.list <- Rarefaction(ind.data,
                                 cov,
                                 function(x, y) RandomSkewers(x, y, iterations)[1],
                                 num.reps=num.reps,
                                 num.cores = num.cores)
  return(comparison.list)
}

#' @rdname Rarefaction
#' @export
RarefactionMantelCor <- function(ind.data,
                                 iterations = 1,
                                 num.reps = 10,
                                 num.cores = 1){
  comparison.list <- Rarefaction(ind.data,
                                 cor,
                                 function(x, y) MantelCor(x, y, iterations)[1],
                                 num.reps=num.reps,
                                 num.cores = num.cores)
  return(comparison.list)
}

#' @rdname Rarefaction
#' @export
RarefactionKrzCor <- function(ind.data,
                              correlation = FALSE,
                              ret.dim = NULL,
                              num.reps = 10,
                              num.cores = 1){
  if(correlation) StatFunc <- cor
  else StatFunc <- cov
  comparison.list <- Rarefaction(ind.data,
                                 StatFunc,
                                 function(x, y) KrzCor(x, y, ret.dim),
                                 num.reps=num.reps,
                                 num.cores = num.cores)
  return(comparison.list)
}
