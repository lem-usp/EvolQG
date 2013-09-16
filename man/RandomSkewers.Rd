\name{RandomSkewers}
\alias{RandomSkewers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RandomSkewers(cov.matrix.1, cov.matrix.2, nsk = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.matrix.1}{
%%     ~~Describe \code{cov.matrix.1} here~~
}
  \item{cov.matrix.2}{
%%     ~~Describe \code{cov.matrix.2} here~~
}
  \item{nsk}{
%%     ~~Describe \code{nsk} here~~
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
function (cov.matrix.1, cov.matrix.2, nsk = 10000) 
{
    traits <- dim(cov.matrix.1)[1]
    base.vector <- Normalize(rnorm(traits))
    random.vectors <- array(rnorm(nsk * traits, mean = 0, sd = 1), 
        c(traits, nsk))
    random.vectors <- apply(random.vectors, 2, Normalize)
    dist <- base.vector \%*\% random.vectors
    dz1 <- apply(cov.matrix.1 \%*\% random.vectors, 2, Normalize)
    dz2 <- apply(cov.matrix.2 \%*\% random.vectors, 2, Normalize)
    real <- apply(dz1 * dz2, 2, sum)
    ac <- mean(real)
    stdev <- sd(real)
    prob <- sum(ac < dist)/nsk
    output <- c(ac, prob, stdev)
    names(output) <- c("AC", "P", "SD")
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
