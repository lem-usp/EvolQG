\name{KrzCor}
\alias{KrzCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Krzanowsky Correlation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KrzCor(cov.matrix.1, cov.matrix.2, ret.dim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.matrix.1}{
%%     ~~Describe \code{cov.matrix.1} here~~
}
  \item{cov.matrix.2}{
%%     ~~Describe \code{cov.matrix.2} here~~
}
  \item{ret.dim}{
%%     ~~Describe \code{ret.dim} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cov.matrix.1, cov.matrix.2, ret.dim = NULL) 
{
    if (is.null(ret.dim)) 
        ret.dim = dim(cov.matrix.1)[1]/2 - 1
    EigenVectors <- function(x) return(eigen(x)$vectors[, 1:ret.dim])
    A <- EigenVectors(cov.matrix.1)
    B <- EigenVectors(cov.matrix.2)
    S <- t(A) \%*\% B \%*\% t(B) \%*\% A
    SL <- sum(eigen(S)$values)/ret.dim
    return(SL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
