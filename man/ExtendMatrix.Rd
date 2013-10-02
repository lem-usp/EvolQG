\name{ExtendMatrix}
\alias{ExtendMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Control Inverse matrix noise with Extension
}
\description{
Calculates the extented covariance matrix estimation
}
\usage{
ExtendMatrix(cov.matrix, cut.off = NULL)
}
\arguments{
  \item{cov.matrix}{
Covariance matrix
}
  \item{cut.off}{
Cut off for second derivative variance. When below first below cut off, last eigen value is selected.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Extended covariance matrix.
}
\references{
Marroig, G., Melo, D. A. R., and Garcia, G. (2012). Modularity, noise, and natural selection. Evolution; international journal of organic evolution, 66(5), 1506-24. doi:10.1111/j.1558-5646.2011.01555.x
}
\author{
Diogo Melo
}
\note{
Covariance matrixe being extended must be larger then 10x10
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
      cov.matrix = RandomMatrix(11, 1, 1, 100)
      ext.matrix = ExtendMatrix(cov.matrix, cut.off = 10e-4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{extension}
\keyword{covariancematrix}
