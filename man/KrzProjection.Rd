\name{KrzProjection}
\alias{KrzProjection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Krzanowsky modified projections
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KrzProjection(cov.matrix.1, cov.matrix.2, ret.dim.1 = NULL, ret.dim.2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.matrix.1}{
%%     ~~Describe \code{cov.matrix.1} here~~
}
  \item{cov.matrix.2}{
%%     ~~Describe \code{cov.matrix.2} here~~
}
  \item{ret.dim.1}{
%%     ~~Describe \code{ret.dim.1} here~~
}
  \item{ret.dim.2}{
%%     ~~Describe \code{ret.dim.2} here~~
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
function (cov.matrix.1, cov.matrix.2, ret.dim.1 = NULL, ret.dim.2 = NULL) 
{
    require(plyr)
    num.traits <- dim(cov.matrix.1)[1]
    if (is.null(ret.dim.1)) 
        ret.dim.1 <- num.traits/2 - 1
    if (is.null(ret.dim.2)) 
        ret.dim.2 <- num.traits/2 - 1
    eigen.cov.matrix.1 <- eigen(cov.matrix.1)
    eVal.1 <- eigen.cov.matrix.1$values
    eVec.1 <- eigen.cov.matrix.1$vectors
    eVar.1 <- t(aaply(1:num.traits, 1, function(n) eVec.1[, n] * 
        sqrt(eVal.1[n])))
    eVec.2 <- eigen(cov.matrix.2)$vectors
    SumSq <- function(x) sum(x^2)
    MapProjection <- function(x) SumSq(t(aaply(1:ret.dim.2, 1, 
        function(n) eVar.1[, x] \%*\% eVec.2[, n])))
    ProjectionNorms <- aaply(1:ret.dim.1, 1, MapProjection)
    output <- list(total.variation = sum(ProjectionNorms)/sum(eVal.1), 
        per.PC = ProjectionNorms/eVal.1[1:ret.dim.1])
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
