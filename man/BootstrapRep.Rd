\name{BootstrapRep}
\alias{BootstrapRep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BootstrapRep(ind.data, nb = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ind.data}{
%%     ~~Describe \code{ind.data} here~~
}
  \item{nb}{
%%     ~~Describe \code{nb} here~~
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
function (ind.data, nb = 100) 
{
    n.ind <- dim(ind.data)[1]
    original.cov.matrix <- var(ind.data)
    v.rep <- c()
    for (N in 1:nb) {
        sampled.data <- sample(1:n.ind, n.ind, TRUE)
        sampled.data.cov.matrix <- var(ind.data[sampled.data, 
            ])
        v.rep[N] <- RandomSkewers(original.cov.matrix, sampled.data.cov.matrix, 
            1000)[1]
    }
    out <- mean(v.rep)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
