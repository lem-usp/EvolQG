BootstrapRep <- function (ind.data, iterations,
                          ComparisonFunc, StatFunc,
                          num.cores)
  # Calculates the repeatability of the covariance matrix of the suplied data
  # via bootstrap ressampling
  #
  # Args:
  #     ind.data: original individual data
  #     nb = number of resamples
  # Return:
  #     returns the mean repeatability
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
  sample.size <-  dim (ind.data) [1]
  c.matrix <- StatFunc(ind.data)
  populations  <- alply(1:iterations, 1,
                        function(x) ind.data[sample (1:sample.size, sample.size, TRUE),],
                        .parallel = parallel)
  comparisons <- laply (populations, function (x) ComparisonFunc (c.matrix, StatFunc(x)),
                        .parallel = parallel)
  return (mean(comparisons))
}

BootstrapRepRandomSkewers <- function(ind.data, iterations = 1000, num.cores = 1){
  repeatability <- BootstrapRep(ind.data, iterations,
                                ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
                                StatFunc = cov,
                                num.cores = num.cores)
  return(repeatability)
}

BootstrapRepMantelCor <- function(ind.data, iterations = 1000, num.cores = 1){
  repeatability <- BootstrapRep(ind.data, iterations,
                                ComparisonFunc = function(x, y) MantelCor(x, y, 1)[1],
                                StatFunc = cor,
                                num.cores = num.cores)
  return(repeatability)
}

BootstrapRepKrzCor <- function(ind.data, iterations = 1000, correlation = F, num.cores = 1){
  if(correlation)  StatFunc <- cor
  else             StatFunc <- cov
  repeatability <- BootstrapRep(ind.data, iterations,
                                ComparisonFunc = function(x, y) MantelCor(x, y, 1)[1],
                                StatFunc = StatFunc,
                                num.cores = num.cores)
  return(repeatability)
}
