\name{MantelCor}
\alias{MantelCor}
\alias{MantelCor.default}
\alias{MantelCor.list}

\title{
    Compare matrices via Mantel Correlation
}
\description{
    Calculates correlation matrix correlation via Mantel.
}
\usage{
MantelCor(cor.matrix, ...)
### Default S3 method
\method{MantelCor}{default}(cor.matrix.1, cor.matrix.2, iterations = 100, mod = FALSE)
### List S3 method
\method{MantelCor}{list}(matrix.list, cor.matrix.2 = NULL,
    repeat.vector = NULL, iterations = 100, mod = FALSE,
    num.cores = 1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cor.matrix, cor.matrix.1, matrix.list}{
        Single covariance matrix or list of covariance matrices.

        If single matrix is suplied, it is compared to cor.matrix.2.

        If list is suplied and no cor.matrix.2 is suplied, all matrices
        are compared.

        If cor.matrix.2 is suplied, all matrices in list are compared to it.
    }
    \item{...}{aditional arguments passed to other methods}
    \item{cor.matrix.2}{
        First argument is compared to cor.matrix.2.
        Optional if cor.matrix.1 is a list.
    }
    \item{iterations}{
        Number of permutations used in significance calculation.
    }
    \item{mod}{
        Set TRUE to use mantel in testing modularity hipotesis. Will return
        AVG+, AVG- and AVG Ratio based on binary hipotesis matrix.
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
    If cor.matrix.1 and cor.matrix.2 are passed, returns matrix pearson
    correlation and significance via mantel permutations.

    If matrix.list and cor.matrix.2 are passed, same as above, but for all matrices in matrix.list.

    If only matrix.list is passed, a matrix of MantelCor average
    values and probabilities of all comparisons.
    If repeat.vector is passed, comparison matrix is corrected below
    diagonal and repeatabilities returned in diagonal.
}
\references{
    http://en.wikipedia.org/wiki/Mantel_test
}
\author{
    Diogo Melo, Guilherme Garcia
}
\note{
    %%  ~~further notes~~
}

\seealso{
    \code{\link{KrzCor}},\code{\link{RandomSkewers}},\code{\link{mantel}},\code{\link{RandomSkewers}},\code{\link{TestModularity}}

}
\examples{
c1 <- RandomMatrix(10)
c2 <- RandomMatrix(10)
c3 <- RandomMatrix(10)
MantelCor(c1, c2)

MantelCor(list(c1, c2, c3))

reps <- unlist(lapply(list(c1, c2, c3), MonteCarloRepMantelCor, 10, 10))
MantelCor(list(c1, c2, c3), repeat.vector = reps)

c4 <- RandomMatrix(10)
MantelCor(list(c1, c2, c3), c4)
}
\keyword{matrixcomparison}
\keyword{matrixcorrelation}
\keyword{randomskewers}
