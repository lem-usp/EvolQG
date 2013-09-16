\name{NotSoRandomMatrix}
\alias{NotSoRandomMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NotSoRandomMatrix(size, first.evec = NULL, eigen.values = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
%%     ~~Describe \code{size} here~~
}
  \item{first.evec}{
%%     ~~Describe \code{first.evec} here~~
}
  \item{eigen.values}{
%%     ~~Describe \code{eigen.values} here~~
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
function (size, first.evec = NULL, eigen.values = NULL) 
{
    ei.vecs <- array(rnorm(size^2), c(size, size))
    if (!is.null(first.evec)) 
        ei.vecs[, 1] = first.evec
    ei.vecs <- apply(ei.vecs, 2, Normalize)
    if (!is.null(eigen.values)) 
        ei.vals <- sort(eigen.values, TRUE)
    else {
        ei.vals <- sort(rchisq(size, 1, 0), TRUE)
        ei.vals <- ei.vals/mean(ei.vals)
    }
    ei.vecs <- qr.Q(qr(ei.vecs))
    rand.mat <- ei.vecs \%*\% diag(ei.vals) \%*\% t(ei.vecs)
    return(rand.mat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
