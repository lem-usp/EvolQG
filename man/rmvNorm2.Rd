\name{rmvNorm2}
\alias{rmvNorm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rmvNorm2(n, theta = rep(0, nrow(sigma)), sigma = diag(length(theta)), method = c("eigen", "svd", "chol"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (n, theta = rep(0, nrow(sigma)), sigma = diag(length(theta)), 
    method = c("eigen", "svd", "chol")) 
{
    if (length(theta) != nrow(sigma)) {
        stop("theta and sigma have non-conforming size")
    }
    sigma.aux <- sigma
    dimnames(sigma.aux) <- NULL
    if (!isTRUE(all.equal(sigma.aux, t(sigma.aux)))) {
        warning("sigma is numerically not symmetric")
    }
    method <- match.arg(method)
    if (method == "eigen") {
        sigma.eigen <- eigen(sigma, symmetric = TRUE)
        if (!all(sigma.eigen$values >= -sqrt(.Machine$double.eps) * 
            abs(sigma.eigen$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        random.deviates <- sigma.eigen$vectors \%*\% diag(sqrt(sigma.eigen$values), 
            length(sigma.eigen$values)) \%*\% t(sigma.eigen$vectors)
    }
    else if (method == "svd") {
        sigma.svd <- svd(sigma)
        if (!all(sigma.svd$d >= -sqrt(.Machine$double.eps) * 
            abs(sigma.svd$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        random.deviates <- t(sigma.svd$v \%*\% (t(sigma.svd$u) * 
            sqrt(sigma.svd$d)))
    }
    else if (method == "chol") {
        sigma.chol <- chol(sigma)
        random.deviates <- matrix(rnorm(n * ncol(sigma)), nrow = n) \%*\% 
            sigma.chol
    }
    random.deviates <- sweep(random.deviates, 2, theta, "+")
    return(random.deviates)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
