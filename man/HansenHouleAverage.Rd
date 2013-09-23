\name{HansenHouleAverage}
\alias{HansenHouleAverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Hansen House Statistics
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HansenHouleAverage(mat, nsk = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
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
function (mat, nsk = 10000) 
{
    with(load.of.functions, {
        n.char <- dim(mat)[1]
        beta.mat <- array(rnorm(n.char * nsk), c(n.char, nsk))
        beta.mat <- apply(beta.mat, 2, Normalize)
        iso.vec <- Normalize(rep(1, times = n.char))
        null.dist <- abs(t(iso.vec) \%*\% beta.mat)
        null.dist <- sort(null.dist)
        crit.value <- null.dist[round(0.95 * nsk)]
        cat("critical value: ", crit.value, "\n")
        parm.dist <- array(0, c(nsk, 8))
        HansenHouleWrap <- function(hh.func) return(apply(beta.mat, 
            2, hh.func, cov.matrix = mat))
        parm.dist[, 1:6] <- sapply(hansen.houle, HansenHouleWrap)
        parm.dist[, 7] <- as.numeric(parm.dist[, 5] > crit.value)
        parm.dist[, 8] <- as.numeric(parm.dist[, 6] > crit.value)
        parm.dist <- cbind(parm.dist, null.dist)
        colnames(parm.dist) <- c("resp", "evol", "cond.evol", 
            "auto", "flex", "const", "flex.n", "const.n", "null.dist")
        parm.av <- colMeans(parm.dist)
        parm.av[7:8] <- parm.av[7:8] * nsk
        pc1 <- eigen(mat)$vectors[, 1]
        HansenHouleWrapPc1 <- function(hh.func) return(hh.func(beta = pc1, 
            cov.matrix = mat))
        maximum <- sapply(hansen.houle, HansenHouleWrapPc1)
        integration <- c(MeanSquaredCorrelation(mat), Pc1Percent(mat))
        names(integration) <- c("MeanSquaredCorrelation", "pc1\%")
        parm.av <- c(integration, parm.av)
        return(list(dist = parm.dist, mean = parm.av, max.val = maximum))
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
