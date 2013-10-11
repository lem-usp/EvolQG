\name{MonteCarloStat}
\alias{MonteCarloR2}
\alias{MonteCarloRepKrzCor}
\alias{MonteCarloRepMantelCor}
\alias{MonteCarloRepRandomSkewers}
\alias{MonteCarloStat}
\title{Parametric repeatabilities with covariance or correlation matrices}
\usage{
  MonteCarloStat(cov.matrix, sample.size, iterations,
    ComparisonFunc, StatFunc, num.cores = 1)

  MonteCarloRepRandomSkewers(cov.matrix, sample.size,
    iterations = 1000, num.cores = 1)

  MonteCarloRepMantelCor(cov.matrix, sample.size,
    iterations = 1000, num.cores = 1)

  MonteCarloRepKrzCor(cov.matrix, sample.size,
    correlation = F, iterations = 1000, num.cores = 1)

  MonteCarloR2(cov.matrix, sample.size, iterations = 1000,
    num.cores = 1)
}
\arguments{
  \item{cov.matrix}{Covariance matrix.}

  \item{sample.size}{Size of the random populations}

  \item{iterations}{Number of random populations}

  \item{ComparisonFunc}{Comparison functions for the
  calculated statistic}

  \item{StatFunc}{Function for calculating the statistic}

  \item{correlation}{When using BootstrapRepKrzCor,
  statistic can be correlation or covariance. If TRUE, uses
  correlation.}

  \item{num.cores}{Number of threads to use in computation.
  Requires doMC library.}
}
\value{
  returns the mean repeatability, or mean value of
  comparisons from samples to original statistic.
}
\description{
  Using a multivariate normal model, random populations are
  generated using the suplied covariance matrix. A
  statistic is calculated on the random population and
  compared to the statistic calculated on the original
  matrix. Prepackaged functions for common comparison
  functions and statistics are suplied.
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


MonteCarloRepRandomSkewers(cov.matrix, 30)
MonteCarloRepMantelCor(cov.matrix, 30)
MonteCarloRepKrzCor(cov.matrix, 30)
MonteCarloRepKrzCor(cov.matrix, 30, TRUE)

#Creating repeatability vector for a list of matrices
mat.list <- RandomMatrix(10, 3, 1, 10)
unlist(lapply(mat.list, MonteCarloRepKrzCor, 30, TRUE))

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

