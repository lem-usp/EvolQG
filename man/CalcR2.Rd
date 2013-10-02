\name{CalcR2}
\alias{CalcR2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Squared Correlations
}
\description{
Calculates the mean squared correlation of a covariance or correlation matrix. Measures integration.
}
\usage{
CalcR2(c.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c.matrix}{
Covariance or correlation matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Mean squared value of of diagonal elements of correlations
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Diogo Melo, Guilherme Garcia
}
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{MeanMatrixStatistics}}, \code{\link{Flexibility}}
}
\examples{
cov.matrix <- RandomMatrix(10, 1, 1, 10)
CalcR2(cov.matrix)
CalcR2(cov2cor(cov.matrix))
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{correlations}
\keyword{integrations}
