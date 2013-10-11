\name{Rarefaction}
\alias{PlotRarefaction}
\alias{Rarefaction}
\alias{RarefactionKrzCor}
\alias{RarefactionMantelCor}
\alias{RarefactionRandomSkewers}
\title{Rarefaction analysis via ressampling
Calculates the repeatability of a statistic of the data, such as
correlation or covariance matrix, via bootstrap ressampling with
varying sample sizes, from 2 to the size of the original data.}
\usage{
  Rarefaction(ind.data, StatFunc, ComparisonFunc,
    num.reps = 10, num.cores = 1)

  PlotRarefaction(comparison.list, y.axis = "Statistic")

  RarefactionRandomSkewers(ind.data, iterations = 1000,
    num.reps = 10, num.cores = 1)

  RarefactionMantelCor(ind.data, iterations = 1,
    num.reps = 10, num.cores = 1)

  RarefactionKrzCor(ind.data, correlation = FALSE,
    ret.dim = NULL, num.reps = 10, num.cores = 1)
}
\arguments{
  \item{ind.data}{Matrix of residuals or indiviual
  measurments}

  \item{num.reps}{number of populations sampled per sample
  size}

  \item{iterations}{Parameter for comparison function.
  Number of random skewers or number of permutations in
  mantel}

  \item{ComparisonFunc}{Comparison function for calculated
  statistic}

  \item{StatFunc}{Statistic to be calculated}

  \item{num.cores}{Number of threads to use in computation.
  Requires doMC library.}

  \item{correlation}{When using RarefactionKrzCor,
  statistic can be correlation or covariance. If TRUE, uses
  correlation.}

  \item{ret.dim}{When using RarefactionKrzCor, number o
  retained dimensions may be specified}

  \item{comparison.list}{output from rarefaction functions
  can be used in ploting}

  \item{y.axis}{Y axis lable in plot}
}
\value{
  returns the mean value of comparisons from samples to
  original statistic, for all sample sizes.
}
\description{
  Samples of various sizes, with replacement, are taken
  from the full population, a statistic calculated and
  compared to the full population statistic. Prepackaged
  functions for common comparison functions and statistics
  are suplied.
}
\details{
  A specialized ploting function displays the results in
  publication quality.

  Bootstraping may be misleading with very small sample
  sizes. Use with caution.
}
\examples{
ind.data <- matrix(rnorm(30*10), 30, 10)
comparison.list <- Rarefaction(ind.data,
                               StatFunc = cov,
                               ComparisonFunc = function(x, y) RandomSkewers(x, y, 1000)[1],
                               num.reps=5,
                               num.cores = 1)

results.RS <- RarefactionRandomSkewers(ind.data, num.reps = 5)
results.Mantel <- RarefactionMantelCor(ind.data, num.reps = 5)
results.KrzCov <- RarefactionKrzCor(ind.data, num.reps = 5)
results.KrzCor <- RarefactionKrzCor(ind.data, TRUE, num.reps = 5)

#Easy access
library(reshape2)
melt(results.RS)

#Plotting using ggplot2
a <- PlotRarefaction(results.RS, "Random Skewers")
b <- PlotRarefaction(results.Mantel, "Mantel")
c <- PlotRarefaction(results.KrzCov, "Krz - Covariance")
d <- PlotRarefaction(results.KrzCor, "Krz - Correlation")

library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(a, vp = vplayout(1, 1))
print(b, vp = vplayout(1, 2))
print(c, vp = vplayout(2, 1))
print(d, vp = vplayout(2, 2))
}
\author{
  Diogo Melo, Guilherme Garcia
}
\seealso{
  \code{\link{BootstrapRep}}
}
\keyword{bootstap}
\keyword{rarefaction}
\keyword{repeatability}

