Rarefaction <- function(ind.data,
                        StatFunc,
                        ComparisonFunc,
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

RarefactionRandomSkewers <- function(ind.data, iterations = 1000, num.reps = 10){
  comparison.list <- Rarefaction(ind.data,
                                 cov,
                                 function(x, y) RandomSkewers(x, y, iterations)[1],
                                 num.rep=num.reps)
  return(comparison.list)
}

RarefactionMantelCor <- function(ind.data, iterations = 1, num.reps = 10){
  comparison.list <- Rarefaction(ind.data,
                                 cor,
                                 function(x, y) MantelCor(x, y, iterations)[1],
                                 num.rep=num.reps)
  return(comparison.list)
}

RarefactionKrzCor <- function(ind.data, correlation = F, ret.dim = NULL, num.reps = 10){
  if(correlation) StatFunc <- cor
  else StatFunc <- cov
  comparison.list <- Rarefaction(ind.data,
                                 StatFunc,
                                 function(x, y) KrzCor(x, y, ret.dim),
                                 num.rep=num.reps)
  return(comparison.list)
}
