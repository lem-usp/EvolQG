\name{TransferSize}
\alias{TransferSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TransferSize(cov.matrix.1, cov.matrix.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.matrix.1}{
%%     ~~Describe \code{cov.matrix.1} here~~
}
  \item{cov.matrix.2}{
%%     ~~Describe \code{cov.matrix.2} here~~
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
function (cov.matrix.1, cov.matrix.2) 
{
    eigen.cov.matrix.1 = eigen(cov.matrix.1)
    eVal.1 = eigen.cov.matrix.1$values
    eVec.1 = eigen.cov.matrix.1$vectors
    eigen.cov.matrix.2 = eigen(cov.matrix.2)
    eVal.2 = eigen.cov.matrix.2$values
    eVec.2 = eigen.cov.matrix.2$vectors
    aux = eVal.1[1]
    eVal.1[1] = eVal.2[1]
    eVal.2[1] = aux
    cov.matrix.1.size.transfer <- eVec.1 \%*\% diag(eVal.1) \%*\% 
        t(eVec.1)
    cov.matrix.2.size.transfer <- eVec.2 \%*\% diag(eVal.2) \%*\% 
        t(eVec.2)
    return(list(cov.matrix.1.size.transfer, cov.matrix.2.size.transfer))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
