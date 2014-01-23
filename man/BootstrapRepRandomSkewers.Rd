\name{BootstrapRepRandomSkewers}
\alias{BootstrapRepRandomSkewers}
\title{Bootstrap analysis via ressampling}
\usage{
BootstrapRepRandomSkewers(ind.data, iterations = 1000, correlation = F,
  num.cores = 1)
}
\arguments{
  \item{ind.data}{Matrix of residuals or indiviual
  measurments}

  \item{iterations}{Number of ressamples to take}

  \item{correlation}{If TRUE, correlation matrix is used,
  else covariance matrix}

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
population, the covariance matrix is calculated and
compared to the full population matrix using Random
Skewers. Prepackaged functions for common comparison
functions and statistics are suplied.
}
\examples{
BootstrapRepRandomSkewers(iris[,1:4], 50)
}
\author{
Diogo Melo, Guilherme Garcia
}
\seealso{
\code{\link{MonteCarloStat}}, \code{\link{AlphaRep}},
\code{\link{BootstrapRepKrzCor}},
\code{\link{BootstrapRepMantelCor}},
\code{\link{BootstrapRep}}
}
\keyword{bootstrap}
\keyword{repetabilities}

