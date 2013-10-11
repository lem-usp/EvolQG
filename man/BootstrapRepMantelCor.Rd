\name{BootstrapRepMantelCor}
\alias{BootstrapRepMantelCor}
\title{Bootstrap analysis via ressampling}
\usage{
  BootstrapRepMantelCor(ind.data, iterations = 1000,
    num.cores = 1)
}
\arguments{
  \item{ind.data}{Matrix of residuals or indiviual
  measurments}

  \item{iterations}{Number of ressamples to take}

  \item{num.cores}{Number of threads to use in computation.
  Requires doMC library.}
}
\value{
  returns the mean repeatability, or mean value of
  comparisons from samples to original statistic.
}
\description{
  Calculates the repeatability of the correlation matrix of
  the suplied data via bootstrap ressampling
}
\details{
  Samples with replacement are taken from the full
  population, the correlation matrix is calculated and
  compared to the full population matrix via Mantel
  correlation
}
\examples{
BootstrapRepMantelCor(iris[,1:4], 50)


#Creating repeatability vector for a list of matrices

mat.list <- RandomMatrix(10, 3, 1, 10)
unlist(lapply(mat.list, BootstrapRepMantelCor), 30, TRUE)
}
\author{
  Diogo Melo, Guilherme Garcia
}
\seealso{
  \code{\link{MonteCarloStat}}, \code{\link{AlphaRep}},
  \code{\link{BootstrapRepKrzCor}},
  \code{\link{BootstrapRep}},
  \code{\link{BootstrapRepRandomSkewers}}
}

