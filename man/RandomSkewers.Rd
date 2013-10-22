\name{RandomSkewers}
\alias{RandomSkewers}
\alias{RandomSkewers.default}
\alias{RandomSkewers.list}
\title{Compare matrices via RandomSkewers}
\usage{
  RandomSkewers(cov.x, cov.y, ...)

  \method{RandomSkewers}{default} (cov.x, cov.y,
    iterations = 1000, ...)

  \method{RandomSkewers}{list} (cov.x, cov.y = NULL,
    iterations = 1000, repeat.vector = NULL, num.cores = 1,
    ...)
}
\arguments{
  \item{cov.x}{Single covariance matrix or list of
  covariance matrices. If single matrix is suplied, it is
  compared to cov.y. If list is suplied and no cov.y is
  suplied, all matrices are compared. If cov.y is suplied,
  all matrices in list are compared to it.}

  \item{cov.y}{First argument is compared to cov.y.
  Optional if cov.x is a list.}

  \item{iterations}{Number of random vectors used in
  comparison.}

  \item{repeat.vector}{Vector of repeatabilities for
  correlation correction.}

  \item{num.cores}{If list is passed, number of threads to
  use in computation. Requires doMC library.}

  \item{...}{aditional arguments passed to other methods}
}
\value{
  If cov.x and cov.y are passed, returns average value of
  response vectors correlation ('correlation'),
  significance ('probability') and standard deviation of
  response vectors correlation ('correlation_sd')

  If cov.x and cov.y are passed, same as above, but for all
  matrices in cov.x.

  If only a list is passed to cov.x, a matrix of
  RandomSkewers average values and probabilities of all
  comparisons. If repeat.vector is passed, comparison
  matrix is corrected above diagonal and repeatabilities
  returned in diagonal.
}
\description{
  Calculates covariance matrix correlation via random
  skewers
}
\examples{
c1 <- RandomMatrix(10)
c2 <- RandomMatrix(10)
c3 <- RandomMatrix(10)
RandomSkewers(c1, c2)

RandomSkewers(list(c1, c2, c3))

reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRepRandomSkewers, 10, 10))
RandomSkewers(list(c1, c2, c3), repeat.vector = reps)

c4 <- RandomMatrix(10)
RandomSkewers(list(c1, c2, c3), c4)
}
\author{
  Diogo Melo, Guilherme Garcia
}
\references{
  Cheverud, J. M., and Marroig, G. (2007). Comparing
  covariance matrices: Random skewers method compared to
  the common principal components model. Genetics and
  Molecular Biology, 30, 461-469.
}
\seealso{
  \code{\link{KrzCor}},\code{\link{MantelCor}}
}
\keyword{matrixcomparison}
\keyword{matrixcorrelation}
\keyword{randomskewers}

