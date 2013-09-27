\name{KrzCor}
\alias{KrzCor}
\alias{KrzCor.default}
\alias{KrzCor.list}

\title{
    Compare matrices via Krzanowski Correlation
}
\description{
    Calculates covariance matrix correlation via random skewers
}
\usage{
KrzCor(cov.matrix, ...)
### Default S3 method
\method{KrzCor}{default}(cov.matrix.1, cov.matrix.2, ret.dim = NULL)

### List S3 method
\method{KrzCor}{list}(matrix.list, cov.matrix.2 = NULL,
       repeat.vector = NULL, ret.dim = NULL,
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
    \item{ret.dim}{
        number of retained dimensions in the comparison,
        default for nxn matrix is n/2-1
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
    Kzranowski correlation
}
\references{
Krzanowski, W. J. (1979). Between-Groups Comparison of Principal
Components. Journal of the American Statistical Association, 74(367),
703. doi:10.2307/2286995
}
\author{
    Diogo Melo, Guilherme Garcia
}
\note{
    %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{RandomSkewers}},\code{\link{KrzProjection}},\code{\link{MantelCor}}
}
\examples{
c1 <- RandomMatrix(10)
c2 <- RandomMatrix(10)
c3 <- RandomMatrix(10)
KrzCor(c1, c2)

KrzCor(list(c1, c2, c3))

reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRepKrzCor, 10, 10))
KrzCor(list(c1, c2, c3), repeat.vector = reps)

c4 <- RandomMatrix(10)
KrzCor(list(c1, c2, c3), c4)
}
\keyword{matrixcomparison}
\keyword{matrixcorrelation}
\keyword{Krzanowski}
