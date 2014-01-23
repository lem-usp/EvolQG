\name{BootstrapRep}
\alias{BootstrapRep}
\title{Bootstrap analysis via ressampling}
\usage{
BootstrapRep(ind.data, iterations, ComparisonFunc, StatFunc, num.cores)
}
\arguments{
  \item{ind.data}{Matrix of residuals or indiviual
  measurments}

  \item{iterations}{Number of ressamples to take}

  \item{ComparisonFunc}{Comparison function for calculated
  statistic}

  \item{StatFunc}{Statistic to be calculated}

  \item{num.cores}{Number of threads to use in computation.
  Requires doMC library.}
}
\value{
returns the mean repeatability, or mean value of
comparisons from samples to original statistic.
}
\description{
Calculates the repeatability of the covariance matrix of
the suplied data via bootstrap ressampling
}
\details{
Samples with replacement are taken from the full
population, a statistic calculated and compared to the full
population statistic. Prepackaged functions for common
comparison functions and statistics are suplied.
}
\examples{
BootstrapRep(iris[,1:4], iterations = 5,
             ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
             StatFunc = cov,
             num.cores = 1)
}
\author{
Diogo Melo, Guilherme Garcia
}
\seealso{
\code{\link{MonteCarloStat}}, \code{\link{AlphaRep}},
\code{\link{BootstrapRepKrzCor}},
\code{\link{BootstrapRepMantelCor}},
\code{\link{BootstrapRepRandomSkewers}}
}
\keyword{bootstrap}
\keyword{repetabilities}

