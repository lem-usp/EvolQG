\name{MultiMahalanobis}
\alias{MultiMahalanobis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Mahalonabis distance for many vectors
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MultiMahalanobis(mean.list, cov.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mean.list}{
%%     ~~Describe \code{mean.list} here~~
}
  \item{cov.matrix}{
%%     ~~Describe \code{cov.matrix} here~~
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
function (mean.list, cov.matrix) 
{
    n.means <- length(mean.list)
    means.names <- names(mean.list)
    distances <- array(0, c(n.means, n.means))
    for (i in 1:(n.means - 1)) {
        for (j in (i + 1):n.means) {
            cat(i, " ", j, "\n")
            comparing.now <- mahalanobis(mean.list[[i]], mean.list[[j]], 
                cov.matrix)
            distances[i, j] <- distances[j, i] <- comparing.now
        }
    }
    rownames(distances) <- means.names
    colnames(distances) <- means.names
    return(distances)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
