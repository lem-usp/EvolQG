\name{Matrix Statistics}
\alias{MeanMatrixStatistics}
\alias{Autonomy}
\alias{ConditionalEvolvability}
\alias{Constraints}
\alias{Evolvability}
\alias{Flexibility}
\alias{Pc1Percent}
\alias{Respondability}
\title{
Calculate mean values for various matrix statistics
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MeanMatrixStatistics (cov.matrix, iterations = 1000,
                      full.results = F, num.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.matrix}{
A covariance matrix
}
  \item{iterations}{
Number of random vectors to be used in calculating the stochastic statistics
}
\item{full.results}{
If TRUE, full distribution of statistics will be returned.
}
\item{num.cores}{
Number of threads to use in calculations. Requires doMC library.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{dist}{Full distribution of stochastic statistics, only if full.resuts == TRUE}
\item{mean}{Mean value for all statistics}
}
\references{
Hansen, T. F., & Houle, D. (2008). Measuring and comparing evolvability
and constraint in multivariate characters. Journal of evolutionary
biology, 21(5), 1201–19. doi:10.1111/j.1420-9101.2008.01573.x
}
\author{
%%  ~~who you are~~
}
\note{
Diogo Melo, Guilherme Garcia
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

cov.matrix <- cov(iris[,1:4])
MeanMatrixStatistics(cov.matrix)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Autonomy}
\keyword{ConditionalEvolvability}
\keyword{Constraints}
\keyword{Evolvability}
\keyword{Flexibility}
\keyword{Pc1Percent}
\keyword{Respondability}
