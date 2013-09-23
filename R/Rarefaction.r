Rarefaction <- function(ind.data,
                        StatFunc = cov,
                        ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
                        num.reps = 10)
{
  library(plyr)
  observed.stat = StatFunc(ind.data)
  num.ind = dim(ind.data)[1]
  MapStatFunc <- function(n){
    SampleFunction <- function(x){
      local.sample = sample(1:num.ind, n, replace=T)
      return(StatFunc(ind.data[local.sample,]))
    }
    return(alply(1:num.reps, 1, SampleFunction))
  }
  sample.stats = alply(2:num.ind, 1, MapStatFunc, .parallel = F)
  MapComparisonFunc <- function(stat.list){
    laply(stat.list, function(x) ComparisonFunc(x, observed.stat))
  }
  comparison.list = llply(sample.stats, MapComparisonFunc, .parallel = F)
  return(comparison.list)
}

PlotRarefaction <- function(comparison.list, y.axis = "Statistic"){
  library(reshape2)
  library(ggplot2)
  plot.df = melt(comparison.list)
  plot.df = as.data.frame(lapply(plot.df, as.numeric))
  rarefaction.plot = ggplot(plot.df, aes(L1, value, group = L1)) +
  layer(geom = "boxplot") +
  scale_x_continuous("ReSample Size") +
  scale_y_continuous(y.axis) +
  theme_bw()
  return(rarefaction.plot)
}
