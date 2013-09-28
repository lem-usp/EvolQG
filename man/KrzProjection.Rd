\name{KrzProjection}
\alias{KrzProjection}
\alias{KrzProjection.default}
\alias{KrzProjection.list}

\title{
    Compare matrices via Modified Krzanowski Correlation
}
\description{
    Calculates the modified Krzanowski correlation between matrices,
    projecting the variance in each principal  components of the first
    matrix in to the ret.dim.2 components of the second for a list of matrices

}
\usage{
KrzProjection(cov.x, cov.y, ...)
### Default S3 method
\method{KrzProjection}{default}(cov.x, cov.y,
              ret.dim.1 = NULL, ret.dim.2 = NULL, ...)

### List S3 method
\method{KrzProjection}{list}(cov.x, cov.y = NULL, ret.dim.1 = NULL, ret.dim.2 = NULL,
              num.cores = 1, full.results = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cov.x}{
        Single covariance matrix ou list of covariance matrices.
        If single matrix is suplied, it is compared to cov.y.
        If list is suplied and no cov.y is suplied, all matrices
        are compared.
        If cov.y is suplied, all matrices in list are compared to it.
    }
    \item{...}{aditional arguments passed to other methods}
    \item{cov.y}{
        First argument is compared to cov.y.
        Ignored if cov.marix.1 is a list.
    }
    \item{ret.dim.1, ret.dim.2}{
        number of retained dimensions in the comparison,
        default for nxn matrix is n/2-1
    }
    \item{num.cores}{
        If list is passed, number of threads to use in computation. 
        Requires doMC library.
    }
    \item{full.results}{
        if FALSE returns only total variance,
        if TRUE also per PC variance.
    }
}
\details{
    %%  ~~ If necessary, more details than the description above ~~
}
\value{
    Ratio of projected variance to total variance, and ratio of projectecd total in each PC
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
    \code{\link{RandomSkewers}},\code{\link{MantelCor}}
}
\examples{
c1 <- RandomMatrix(10)
c2 <- RandomMatrix(10)
c3 <- RandomMatrix(10)
KrzProjection(c1, c2)

KrzProjection(list(c1, c2, c3))
KrzProjection(list(c1, c2, c3), 5, 4)
KrzProjection(list(c1, c2, c3), 4, 5)

}
\keyword{matrixcomparison}
\keyword{matrixcorrelation}
\keyword{Krzanowski}
