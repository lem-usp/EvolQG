\name{MonteCarloStat}
\alias{MonteCarloR2}
\alias{MonteCarloStat}
\title{Parametric population samples with covariance or correlation matrices}
\usage{
MonteCarloStat(cov.matrix, sample.size, iterations, ComparisonFunc, StatFunc,
  num.cores = 1)

MonteCarloR2(cov.matrix, sample.size, iterations = 1000, num.cores = 1)
}
\arguments{
  \item{cov.matrix}{Covariance matrix.}

  \item{sample.size}{Size of the random populations}

  \item{iterations}{Number of random populations}

  \item{ComparisonFunc}{Comparison functions for the
  calculated statistic}

  \item{StatFunc}{Function for calculating the statistic}

  \item{num.cores}{Number of threads to use in computation.
  Requires doMC library.}
}
\value{
returns the mean repeatability, or mean value of
comparisons from samples to original statistic.
}
\description{
Using a multivariate normal model, random populations are
generated using the suplied covariance matrix. A statistic
is calculated on the random population and compared to the
statistic calculated on the original matrix.
}
\details{
Since this function uses multivariate normal model to
generate populations, only covariance matrices should be
used.
}
\examples{
cov.matrix <- RandomMatrix(10, 1, 1, 10)

MonteCarloStat(cov.matrix, sample.size = 30, iterations = 1000,
               ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
               StatFunc = cov,
               num.cores = 1)


#Calculating R2 confidence intervals
r2.dist <- MonteCarloR2(RandomMatrix(10, 1, 1, 10), 30)
quantile(r2.dist)
}
\author{
Diogo Melo Guilherme Garcia
}
\seealso{
\code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
}
\keyword{montecarlo}
\keyword{parametricsampling}
\keyword{repeatability}

