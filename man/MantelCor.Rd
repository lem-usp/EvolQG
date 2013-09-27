\name{MantelCor}
\alias{MantelCor}

\title{
Compare matrices via Mantel Correlation
}
\description{
Calculates correlation matrix correlation via Mantel.
}
\usage{
MantelCor(x, ...)

### Default S3 method
MantelCor(cor.matrix.1, cor.matrix.2, iterations = 100, mod = F)

### List S3 method
MantelCor(matrix.list, cor.matrix = NULL,
              repeat.vector = NULL, iterations = 100, mod = F,
              num.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cor.matrix.1, matrix.list}{
      Single covariance matrix ou list of covariance matrices.

      If single matrix is suplied, it is compared to cor.matrix.2.

      If list is suplied and no cor.matrix.2 is suplied, all matrices
      are compared.

      If cor.matrix.2 is suplied, all matrices in list are compared to it.
}
  \item{cor.matrix.2}{
      Optional if cor.marix.1 is a list.
      First argument is compared to cor.matrix.2.
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
 \code{\link{KrzCor}},\code{\link{RandomSkewers}},\code{\link{TestModularity}}
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
