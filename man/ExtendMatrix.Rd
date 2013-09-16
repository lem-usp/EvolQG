\name{ExtendMatrix}
\alias{ExtendMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ExtendMatrix(cov.matrix, cutoff = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.matrix}{
%%     ~~Describe \code{cov.matrix} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
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
function (cov.matrix, cutoff = NULL) 
{
    p = dim(cov.matrix)[1]
    if (p < 10) 
        stop("matrix is too small")
    eigen.cov.matrix = eigen(cov.matrix)
    eVal = eigen.cov.matrix$values
    eVec = eigen.cov.matrix$vectors
    if (is.null(cutoff)) {
        grad = array(dim = c(p - 2))
        tr.cov.matrix = sum(eVal)
        for (i in 1:(p - 2)) grad[i] = abs(eVal[i]/tr.cov.matrix - 
            2 * (eVal[i + 1]/tr.cov.matrix) + eVal[i + 2]/tr.cov.matrix)
        var.grad = array(dim = c(p - 6))
        for (i in 1:(p - 6)) {
            var.grad[i] = var(grad[i:(i + 4)])
        }
        length(var.grad[var.grad < 1e-04])
        x11()
        plot(4:(p - 3), var.grad)
        cutoff = floor(locator(1)$x)
    }
    eVal[eVal < eVal[cutoff]] = eVal[cutoff]
    extended.cov.matrix = eVec \%*\% diag(eVal) \%*\% t(eVec)
    colnames(extended.cov.matrix) = colnames(cov.matrix)
    rownames(extended.cov.matrix) = rownames(cov.matrix)
    return(extended.cov.matrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
