\name{MonteCarloRep}
\alias{MonteCarloRep}
\alias{MonteCarloRepRandomSkewers}
\alias{MonteCarloRepMantelCor}
\alias{MonteCarloRepKrzCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Parametric repeatabilities with covariance or correlation matrices
}
\description{
    Using a multivariate normal model, random populations are generated
    using the suplied covariance matrix. A statistic is calculated on the
    random population and compared to the statistic calculated on the
    original matrix. Prepackaged functions for common comparison functions
    and statistics are suplied.

}
\usage{
MonteCarloRep (cov.matrix, sample.size, iterations,
               ComparisonFunc, StatFunc,
               num.cores = 1)
MonteCarloRepRandomSkewers (cov.matrix, sample.size, iterations = 1000, num.cores = 1)
MonteCarloRepMantelCor (cov.matrix, sample.size, iterations = 1000, num.cores = 1)
MonteCarloRepKrzCor (cov.matrix, sample.size, correlation = F, iterations = 1000, num.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cov.matrix}{
Covariance matrix.
    }
    \item{sample.size}{
Size of the random populations
    }
    \item{iterations}{
Number of random populations
    }
    \item{ComparisonFunc}{
Comparison functions for the calculated statistic
    }
    \item{StatFunc}{
Function for calculating the statistic
    }
     \item{correlation}{
        When using BootstrapRepKrzCor, statistic can be correlation or covariance. If TRUE, uses correlation.
    }
     \item{num.cores}{
        Number of threads to use in computation. 
        Requires doMC library.
    }
}
\details{
Since this function uses multivariate normal model to generate
popualtions, only covariance matrices should be used.
}
\value{
     returns the mean repeatability, or mean value of comparisons from
samples to original statistic.
}
\references{
    %% ~put references to the literature/web site here ~
}
\author{
Diogo Melo, Guilherme Garcia
}
\note{
    %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
     \code{\link{BootstrapRep}}, \code{\link{AlphaRep}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
cov.matrix <- RandomMatrix(10, 1, 1, 10)

MonteCarloRep(cov.matrix, sample.size = 30, iterations = 1000,
              ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
              StatFunc = cov,
              num.cores = 1)


MonteCarloRepRandomSkewers(cov.matrix, 30)
MonteCarloRepMantelCor(cov.matrix, 30)
MonteCarloRepKrzCor(cov.matrix, 30)
MonteCarloRepKrzCor(cov.matrix, 30, TRUE)

#Creating repeatability vector for a list of matrices
mat.list <- RandomMatrix(10, 3, 1, 10)
unlist(lapply(mat.list, MonteCarloRepKrzCor, 30, TRUE))
}
\keyword{parametricsampling}
\keyword{montecarlo}
\keyword{repeatability}
