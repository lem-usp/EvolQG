\name{TestModularity}
\alias{TestModularity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test modularity hipotesis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TestModularity(cor.matrix, modularity.hipot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cor.matrix}{
%%     ~~Describe \code{cor.matrix} here~~
}
  \item{modularity.hipot}{
%%     ~~Describe \code{modularity.hipot} here~~
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
function (cor.matrix, modularity.hipot) 
{
    no.hip <- dim(modularity.hipot)[2]
    traits <- dim(modularity.hipot)[1]
    m.hip.array <- array(0, c(traits, traits, no.hip + 1))
    for (N in 1:no.hip) {
        for (L in 1:traits) {
            for (M in 1:traits) {
                m.hip.array[L, M, N] <- ifelse(modularity.hipot[L, 
                  N] & modularity.hipot[M, N], 1, 0)
            }
        }
    }
    m.hip.array[, , no.hip + 1] <- as.integer(as.logical(apply(m.hip.array, 
        c(1, 2), sum)))
    no.hip <- no.hip + 1
    output <- array(0, c(5, no.hip))
    for (N in 1:no.hip) {
        tmp <- MantelCor(cor.matrix, m.hip.array[, , N], mod = TRUE)
        output[, N] <- tmp
    }
    dimnames(output) <- list(names(tmp), c(colnames(modularity.hipot), 
        "Full Integration"))
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
