\name{MultiKrzProjection}
\alias{MultiKrzProjection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MultiKrzProjection(mat.list, ret.dim.1 = NULL, ret.dim.2 = NULL, num.cores = 1, full.results = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat.list}{
%%     ~~Describe \code{mat.list} here~~
}
  \item{ret.dim.1}{
%%     ~~Describe \code{ret.dim.1} here~~
}
  \item{ret.dim.2}{
%%     ~~Describe \code{ret.dim.2} here~~
}
  \item{num.cores}{
%%     ~~Describe \code{num.cores} here~~
}
  \item{full.results}{
%%     ~~Describe \code{full.results} here~~
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
function (mat.list, ret.dim.1 = NULL, ret.dim.2 = NULL, num.cores = 1, 
    full.results = F) 
{
    require(plyr)
    require(reshape2)
    if (num.cores > 1) {
        library(doMC)
        library(foreach)
        registerDoMC(num.cores)
        parallel = TRUE
    }
    else parallel = FALSE
    if (full.results) {
        CompareToNProj <- function(n) llply(mat.list, function(x) {
            KrzProjection(x, mat.list[[n]], ret.dim.1, ret.dim.2)
        }, .parallel = parallel)
    }
    else {
        CompareToNProj <- function(n) llply(mat.list[names(mat.list) != 
            n], function(x) {
            KrzProjection(x, mat.list[[n]], ret.dim.1, ret.dim.2)[1]
        }, .parallel = parallel)
    }
    comparisons.proj <- llply(names(mat.list), CompareToNProj, 
        .progress = "text", .parallel = parallel)
    if (full.results) {
        names(comparisons.proj) = names(mat.list)
        return(comparisons.proj)
    }
    else {
        comparisons.proj <- melt(comparisons.proj)
        comparisons.proj[, 4] = names(mat.list)[(comparisons.proj[, 
            4])]
    }
    return(comparisons.proj[, -2])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
