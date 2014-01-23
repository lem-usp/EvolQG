\name{MeanMatrixStatistics}
\alias{Autonomy}
\alias{ConditionalEvolvability}
\alias{Constraints}
\alias{Evolvability}
\alias{Flexibility}
\alias{MeanMatrixStatistics}
\alias{Pc1Percent}
\alias{Respondability}
\title{Calculate mean values for various matrix statistics}
\usage{
MeanMatrixStatistics(cov.matrix, iterations = 1000, full.results = F,
  num.cores = 1)
}
\arguments{
  \item{cov.matrix}{A covariance matrix}

  \item{iterations}{Number of random vectors to be used in
  calculating the stochastic statistics}

  \item{full.results}{If TRUE, full distribution of
  statistics will be returned.}

  \item{num.cores}{Number of threads to use in
  calculations. Requires doMC library.}
}
\value{
dist Full distribution of stochastic statistics, only if
full.resuts == TRUE

mean Mean value for all statistics
}
\description{
Calculates: Mean Squared Correlation, Autonomy,
ConditionalEvolvability, Constraints, Evolvability,
Flexibility, Pc1Percent, Respondability.
}
\examples{
cov.matrix <- cov(iris[,1:4])
MeanMatrixStatistics(cov.matrix)
}
\author{
Diogo Melo Guilherme Garcia
}
\references{
Hansen, T. F., and Houle, D. (2008). Measuring and
comparing evolvability and constraint in multivariate
characters. Journal of evolutionary biology, 21(5),
1201-19. doi:10.1111/j.1420-9101.2008.01573.x
}
\keyword{Autonomy}
\keyword{ConditionalEvolvability}
\keyword{Constraints}
\keyword{Evolvability}
\keyword{Flexibility}
\keyword{Pc1Percent}
\keyword{Respondability}

