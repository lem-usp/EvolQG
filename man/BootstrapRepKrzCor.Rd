\name{BootstrapRepKrzCor}
\alias{BootstrapRepKrzCor}
\title{Bootstrap analysis via ressampling}
\usage{
  BootstrapRepKrzCor(ind.data, iterations = 1000,
    correlation = F, num.cores = 1)
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
  population, a statistic calculated and compared to the
  full population statistic via Krzanowsky correlation.
}
\examples{
BootstrapRepKrzCor(iris[,1:4], 50)

BootstrapRepKrzCor(iris[,1:4], 50, TRUE)

#Creating repeatability vector for a list of matrices

mat.list <- RandomMatrix(10, 3, 1, 10)
unlist(lapply(mat.list, BootstrapRepKrzCor), 30, TRUE)
}
\author{
  Diogo Melo, Guilherme Garcia
}
\seealso{
  \code{\link{MonteCarloStat}}, \code{\link{AlphaRep}},
  \code{\link{BootstrapRep}},
  \code{\link{BootstrapRepMantelCor}},
  \code{\link{BootstrapRepRandomSkewers}}
}
\keyword{bootstrap}
\keyword{repetabilities}

