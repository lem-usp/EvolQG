\name{ParallelRsMantel}
\alias{ParallelRsMantel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ParallelRsMantel(matrix.list, MatrixCompFunc = RandomSkewers, repeat.vector = NULL, iterations = 1000, n.cores = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix.list}{
%%     ~~Describe \code{matrix.list} here~~
}
  \item{MatrixCompFunc}{
%%     ~~Describe \code{MatrixCompFunc} here~~
}
  \item{repeat.vector}{
%%     ~~Describe \code{repeat.vector} here~~
}
  \item{iterations}{
%%     ~~Describe \code{iterations} here~~
}
  \item{n.cores}{
%%     ~~Describe \code{n.cores} here~~
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
function (matrix.list, MatrixCompFunc = RandomSkewers, repeat.vector = NULL, 
    iterations = 1000, n.cores = 2) 
{
    require(multicore)
    n.matrix <- length(matrix.list)
    matrix.names <- names(matrix.list)
    correlations <- array(0, c(n.matrix, n.matrix))
    rownames(correlations) <- matrix.names
    colnames(correlations) <- matrix.names
    probabilities <- correlations
    index <- which(upper.tri(diag(n.matrix)), arr.ind = TRUE)
    singleComp = function(k) {
        i <- index[k, 1]
        j <- index[k, 2]
        cat(i, " ", j, "\n")
        comparing.now <- MatrixCompFunc(matrix.list[[i]], matrix.list[[j]], 
            iterations)
        return(comparing.now)
    }
    corr.list = mclapply(1:nrow(index), singleComp, mc.cores = n.cores)
    for (k in 1:length(corr.list)) {
        i <- index[k, 1]
        j <- index[k, 2]
        correlations[i, j] <- corr.list[[k]][1]
        probabilities[i, j] <- corr.list[[k]][2]
        if (!is.null(repeat.vector)) 
            correlations[j, i] <- correlations[i, j]/sqrt(repeat.vector[i] * 
                repeat.vector[j])
    }
    if (!is.null(repeat.vector)) {
        diag(correlations) <- repeat.vector
    }
    output <- list(correlations = correlations, probabilities = probabilities)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
