#' Rarefaction analysis via ressampling
#' 
#' Calculates the repeatability of a statistic of the data, such as
#' correlation or covariance matrix, via bootstrap resampling with
#' varying sample sizes, from 2 to the size of the original data.
#'
#' Samples of various sizes, with replacement, are taken from the full population, a statistic calculated
#' and compared to the full population statistic. 
#'
#' A specialized ploting function displays the results in publication quality.
#' @aliases  Rarefaction PlotRarefaction
#' @param ind.data Matrix of residuals or indiviual measurments
#' @param ComparisonFunc comparison function
#' @param ... Aditional arguments passed to ComparisonFunc
#' @param num.reps number of populations sampled per sample size
#' @param iterations Parameter for comparison function. Number of random skewers or number of permutations in mantel.
#' @param num.cores Number of threads to use in computation. The doMC library must be loaded.
#' @param correlation If TRUE, correlation matrix is used, else covariance matrix. MantelCor always uses correlation matrix.
#' @param ret.dim When using Krzanowski Correlation, number of retained dimensions may be specified
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
#' ind.data <- iris[,1:4]
#' 
#' results.RS <- Rarefaction(ind.data, PCAsimilarity, num.reps = 5)
#' results.Mantel <- Rarefaction(ind.data, MatrixCor, correlation = TRUE, num.reps = 5)
#' results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5)
#' results.KrzCor <- Rarefaction(ind.data, KrzCor, correlation = TRUE, num.reps = 5)
#' 
#' #Multiple threads can be used with doMC library
#' library(doMC)
#' results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5, num.cores = 2)
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
                        ComparisonFunc,
                        ...,
                        iterations = 1000, 
                        num.reps = 10,
                        correlation = FALSE, 
                        ret.dim = NULL,
                        num.cores = 1){
  if(correlation)  {StatFunc <- cor; c2v <- cov2cor
  } else {StatFunc <- cov; c2v <- function(x) x}
  rarefaction.list <- Rarefaction_primitive(ind.data,
                                            StatFunc = StatFunc,
                                            ComparisonFunc = function(x, y) ComparisonFunc(c2v(x), 
                                                                                           c2v(y), ...),
                                            num.reps = num.reps,
                                            num.cores = num.cores)
  return(rarefaction.list)
}

Rarefaction_primitive <- function(ind.data,
                                  StatFunc,
                                  ComparisonFunc,
                                  num.reps = 10,
                                  num.cores = 1)
{
  if (num.cores > 1) {
    doMC::registerDoMC(num.cores)
    parallel = TRUE
  } else{
    parallel = FALSE
  }
  if(isSymmetric(as.matrix(ind.data))) stop("input appears to be a matrix, use residuals.")
  observed.stat = StatFunc(ind.data)
  num.ind = dim(ind.data)[1]
  MapStatFunc <- function(n){
    SampleFunction <- function(x){
      while(TRUE){
        local.sample = sample(1:num.ind, n, replace=T)
        out <- tryCatch(StatFunc(ind.data[local.sample,]), warning=function(w) w)
        if(!is(out, "warning"))
          return(out)
      }
    }
    return(alply(1:num.reps, 1, SampleFunction))
  }
  sample.stats = alply(2:num.ind, 1, MapStatFunc, .parallel = parallel)
  MapComparisonFunc <- function(stat.list){
    ldply(stat.list, function(x) ComparisonFunc(x, observed.stat))[,2]
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
  scale_x_continuous("Resample Size") +
  scale_y_continuous(y.axis) +
  theme_bw()
  return(rarefaction.plot)
}
