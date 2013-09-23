\name{MantelCor}
\alias{MantelCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix Correlation with Mantel significance
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MantelCor(cor.matrix.1, cor.matrix.2, nit = 1000, mod = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cor.matrix.1}{
%%     ~~Describe \code{cor.matrix.1} here~~
}
  \item{cor.matrix.2}{
%%     ~~Describe \code{cor.matrix.2} here~~
}
  \item{nit}{
%%     ~~Describe \code{nit} here~~
}
  \item{mod}{
%%     ~~Describe \code{mod} here~~
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
function (cor.matrix.1, cor.matrix.2, nit = 1000, mod = FALSE) 
{
    fixed.matrix <- cor.matrix.1[lower.tri(cor.matrix.1)]
    shuffled.matrix <- cor.matrix.2[lower.tri(cor.matrix.2)]
    index <- shuffled.matrix
    correlation <- cor(fixed.matrix, shuffled.matrix)
    shuffled.correlation <- c()
    for (N in 1:nit) {
        shuffle <- sample(1:dim(cor.matrix.1)[1])
        shuffled.matrix <- cor.matrix.2[shuffle, shuffle][lower.tri(cor.matrix.1)]
        shuffled.correlation[N] <- cor(fixed.matrix, shuffled.matrix)
    }
    prob <- sum(shuffled.correlation > as.vector(correlation))/nit
    if (mod == TRUE) {
        avg.plus <- mean(cor.matrix.1[lower.tri(cor.matrix.1)][index != 
            0])
        avg.minus <- mean(cor.matrix.1[lower.tri(cor.matrix.1)][index == 
            0])
        avg.ratio <- avg.plus/avg.minus
        output <- c(correlation, prob, avg.plus, avg.minus, avg.ratio)
        names(output) <- c("R\302\262", "Probability", "AVG+", 
            "AVG-", "AVG Ratio")
    }
    else {
        output <- c(correlation, prob)
        names(output) <- c("R\302\262", "Probability")
    }
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
