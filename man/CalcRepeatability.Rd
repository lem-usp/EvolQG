\name{CalcRepeatability}
\alias{CalcRepeatability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Repeatabilities
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CalcRepeatability(ID, ind.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
%%     ~~Describe \code{ID} here~~
}
  \item{ind.data}{
%%     ~~Describe \code{ind.data} here~~
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
function (ID, ind.data) 
{
    models.list <- apply(ind.data, 2, function(vec) {
        return(lm(vec ~ ID))
    })
    models.list <- lapply(models.list, anova)
    rep.itself <- function(lm.model) {
        msq <- lm.model$"Mean Sq"
        s2a <- (msq[1] - msq[2])/2
        out <- s2a/(s2a + msq[2])
        return(out)
    }
    out <- sapply(models.list, rep.itself)
    names(out) <- colnames(ind.data)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
