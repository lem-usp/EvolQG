\name{RemoveSize}
\alias{RemoveSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Remove Size Variation
}
\description{
    Removes first principal component effect in a covariance matrix.
}
\usage{
    RemoveSize(cov.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cov.matrix}{
        Covariance matrix
    }
}
\details{
    Function sets the first eigen value to zero.
}
\value{
    Altered covariance matriz, with no variation on former first principal component
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
    cov.matrix <- RandomMatrix(10, 1, 1, 10)
    no.Size.cov.matrix <- RemoveSize(cov.matrix)
    eigen(cov.matrix)
    eigen(no.size.cov.matrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{size}
