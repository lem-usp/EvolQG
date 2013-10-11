\name{RandCorr}
\alias{RandCorr}
\title{Random correlation matrix}
\usage{
  RandCorr(num.traits, ke = 10^-3)
}
\arguments{
  \item{num.traits}{Number of traits in random matrix}

  \item{ke}{Parameter for correlation matrix generation.
  Involves check for positive defitness}
}
\value{
  Random Matrix
}
\description{
  Provides random correlation matrices for quick tests.
  Should not be used for statistics.
}
\examples{
#single 10x10 correlation matrix
RandomCorr(10)
}
\author{
  Diogo Melo Edgar Zanella
}

