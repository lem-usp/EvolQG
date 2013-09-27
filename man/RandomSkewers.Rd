\name{RandomSkewers}
\alias{RandomSkewers}
\alias{RandomSkewers.default}
\alias{RandomSkewers.list}

\title{
    Compare matrices via RandomSkewers
}
\description{
    Calculates covariance matrix correlation via random skewers
}
\usage{
RandomSkewers(cov.matrix, ...)
\method{RandomSkewers}{default}(cov.matrix.1, cov.matrix.2, iterations = 1000)
\method{RandomSkewers}{list}(matrix.list, cov.matrix.2 = NULL,
              repeat.vector = NULL, iterations = 1000, 
              num.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cov.matrix, cov.matrix.1, matrix.list}{
        Single covariance matrix or list of covariance matrices.
        If single matrix is suplied, it is compared to cov.matrix.2.
        If list is suplied and no cov.matrix.2 is suplied, all matrices
        are compared.
        If cov.matrix.2 is suplied, all matrices in list are compared to it.
    }
    \item{...}{aditional arguments passed to other methods}
    \item{cov.matrix.2}{
        First argument is compared to cov.matrix.2.
        Optional if cov.matrix.1 is a list.
    }
    \item{iterations}{
        Number of random vectors used in comparison.
    }
    \item{repeat.vector}{
        Vector of repeatabilities for correlation correction.
    }
    \item{num.cores}{
        If list is passed, number of threads to use in computation. 
        Requires doMC library.
    }
}
\details{
    %%  ~~ If necessary, more details than the description above ~~
}
\value{
    If cov.matrix.1 and cov.matrix.2 are passed, returns average value
    of random skewers ('AC'), significance ('prob') and standard deviation
    of random skewers ('SD')

    If matrix.list and cov.matrix.2 are passed, same as above, but for all matrices in matrix.list.

    If only matrix.list is passed, a matrix of RandomSkewers average
    values and probabilities of all comparisons.
    If repeat.vector is passed, comparison matrix is corrected below
    diagonal and repeatabilities returned in diagonal.
}
\references{
    Cheverud, J. M., and Marroig, G. (2007). Comparing covariance matrices:
    Random skewers method compared to the common principal components model.
    Genetics and Molecular Biology, 30, 461-469.
}
\author{
    Diogo Melo, Guilherme Garcia
}
\note{
    %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{KrzCor}},\code{\link{MantelCor}}
}
\examples{
c1 <- RandomMatrix(10)
c2 <- RandomMatrix(10)
c3 <- RandomMatrix(10)
RandomSkewers(c1, c2)

RandomSkewers(list(c1, c2, c3))

reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRepRandomSkewers, 10, 10))
RandomSkewers(list(c1, c2, c3), repeat.vector = reps)

c4 <- RandomMatrix(10)
RandomSkewers(list(c1, c2, c3), c4)
}
\keyword{matrixcomparison}
\keyword{matrixcorrelation}
\keyword{randomskewers}
