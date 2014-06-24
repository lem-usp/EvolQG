\name{ExtendMatrix}
\alias{ExtendMatrix}
\title{Control Inverse matrix noise with Extension}
\usage{
ExtendMatrix(cov.matrix, cut.off = NULL)
}
\arguments{
  \item{cov.matrix}{Covariance matrix}

  \item{cut.off}{Cut off for second derivative variance.
  When below first below cut off, last eigen value is
  selected.}
}
\value{
Extended covariance matrix.
}
\description{
Calculates the extented covariance matrix estimation
}
\note{
Covariance matrix being extended must be larger then 10x10
}
\examples{
cov.matrix = RandomMatrix(11, 1, 1, 100)
ext.matrix = ExtendMatrix(cov.matrix, cut.off = 10e-4)
}
\author{
Diogo Melo
}
\references{
Marroig, G., Melo, D. A. R., and Garcia, G. (2012).
Modularity, noise, and natural selection. Evolution;
international journal of organic evolution, 66(5), 1506-24.
doi:10.1111/j.1558-5646.2011.01555.x
}
\keyword{covariancematrix}
\keyword{extension}

