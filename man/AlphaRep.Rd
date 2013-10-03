\name{AlphaRep}
\alias{AlphaRep}

\title{
    Alpha Repetability
}
\description{
    Calculates the matrix repeatability using the equation in Cheverud 1996
    Quantitative genetic analysis of cranial morphology in the cotton-top
    (Saguinus oedipus) and saddle-back (S. fuscicollis) tamarins. Journal of Evolutionary Biology 9, 5-42.
}
\usage{
AlphaRep(cor.matrix, sample.size)
}

\arguments{
    \item{cor.matrix}{
        Correlation matrix
    }
    \item{sample.size}{
        Sample size used in matrix estimatio
    }
}
\details{
    %%  ~~ If necessary, more details than the description above ~~
}
\value{
    Alpha repeatability
}
\references{
    Cheverud 1996 Quantitative genetic analysis of cranial morphology in the cotton-top
    (Saguinus oedipus) and saddle-back (S. fuscicollis) tamarins. Journal of Evolutionary Biology 9, 5-42.
}
\author{
    Diogo Melo, Guilherme Garcia
}
\note{
    %%  ~~further notes~~
}

\seealso{
    \code{\link{MonteCarloStat}}, \code{\link{BootstrapRep}}
}
\examples{

#For single matrices
cor.matrix <- RandomMatrix(10)
AlphaRep(cor.matrix, 10)
AlphaRep(cor.matrix, 100)

#For many matrices
mat.list <- RandomMatrix(10, 100)
sample.sizes <- floor(runif(100, 20, 50))
mat.plus.samples <- alply(1:100, 1, function(x) list(mat.list[[x]], sample.sizes[x]))
laply(mat.plus.samples, function(x) AlphaRep(x[[1]], x[[2]]))
}

\keyword{repetatability}
