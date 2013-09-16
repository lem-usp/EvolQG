\name{MultiKrzCor}
\alias{MultiKrzCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MultiKrzCor(matrix.list, repeat.vector = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix.list}{
%%     ~~Describe \code{matrix.list} here~~
}
  \item{repeat.vector}{
%%     ~~Describe \code{repeat.vector} here~~
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
function (matrix.list, repeat.vector = NULL) 
{
    n.matrix <- length(matrix.list)
    matrix.names <- names(matrix.list)
    correlations <- array(0, c(n.matrix, n.matrix))
    for (i in 1:(n.matrix - 1)) {
        for (j in (i + 1):n.matrix) {
            cat(i, " ", j, "\n")
            comparing.now <- KrzCor(matrix.list[[i]], matrix.list[[j]])
            correlations[i, j] <- comparing.now
            if (!is.null(repeat.vector)) 
                correlations[j, i] <- correlations[i, j]/sqrt(repeat.vector[i] * 
                  repeat.vector[j])
        }
    }
    if (!is.null(repeat.vector)) {
        diag(correlations) <- repeat.vector
    }
    rownames(correlations) <- matrix.names
    colnames(correlations) <- matrix.names
    return(correlations)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
