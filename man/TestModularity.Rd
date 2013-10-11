\name{TestModularity}
\alias{TestModularity}
\title{Test modularity hipotesis}
\usage{
  TestModularity(cor.matrix, modularity.hipot,
    iterations = 100)
}
\arguments{
  \item{cor.matrix}{Correlation matrix}

  \item{modularity.hipot}{modularity.hipot: Matrix of
  hipotesis. Each line represents a trait and each column a
  module. if modularity.hipot[i,j] == 1, trait i is in
  module j.}

  \item{iterations}{Number of iterations, to be passed to
  MantelCor}
}
\value{
  Returns mantel correlation and associated probability for
  each modularity hipotesis, along with AVG+, AVG-, AVG
  Ratio for each module. A total hipotesis combining all
  hipotesis is also tested.
}
\description{
  Tests modularity hipotesis using cor.matrix matrix and
  trait groupings
}
\examples{
cor.matrix <- RandomMatrix(10)
rand.hipots <- matrix(sample(c(1, 0), 30, replace=TRUE), 10, 3)
mod.test <- TestModularity(cor.matrix, rand.hipots)
}
\author{
  Diogo Melo, Guilherme Garcia
}
\seealso{
  \code{\link{MantelCor}}
}
\keyword{mantel}
\keyword{modularity}

