\name{MonteCarloR2}
\alias{MonteCarloR2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MonteCarloR2(cov.matrix, sample.size, iterations = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.matrix}{
%%     ~~Describe \code{cov.matrix} here~~
}
  \item{sample.size}{
%%     ~~Describe \code{sample.size} here~~
}
  \item{iterations}{
%%     ~~Describe \code{iterations} here~~
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
function (cov.matrix, sample.size, iterations = 1000) 
{
    n.traits <- dim(cov.matrix)[1]
    populations <- list()
    for (i in 1:iterations) populations[[i]] <- rmvNorm2(sample.size, 
        sigma = cov.matrix, method = "chol")
    it.matrices <- lapply(populations, cor)
    it.r2 <- sapply(it.matrices, CalcR2)
    r2 <- c(CalcR2(cov2cor(cov.matrix)), it.r2)
    return(r2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
