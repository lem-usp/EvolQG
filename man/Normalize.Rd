\name{Normalize}
\alias{Normalize}
\alias{Norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalize and Norm
}
\description{
Norm returns the euclidian norm of a vector, normalze returns a vector with unit norm.
}
\usage{
Normalize(x)
Norm(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    Normalized vector or inpout vector norm.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
    x <- rnorm(10)
    n.x <- Normalize(x)
    Norm(x)
    Norm(n.x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
