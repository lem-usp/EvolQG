\name{BootstrapRep}
\alias{BootstrapRep}
\alias{BootstrapRepRandomSkewers}
\alias{BootstrapRepMantelCor}
\alias{BootstrapRepKrzCor}
\title{
    Calculates the repeatability of the covariance matrix of the suplied data
    via bootstrap ressampling
}
\description{
    Samples with replacement are taken from the full population, a statistic calculated
    and compared to the full population statistic. Prepackaged functions
    for common comparison functions and statistics are suplied.
}
\usage{
BootstrapRep(ind.data, iterations, ComparisonFunc, StatFunc, num.cores)
BootstrapRepRandomSkewers(ind.data, iterations = 1000, num.cores = 1)
BootstrapRepMantelCor(ind.data, iterations = 1000, num.cores = 1)
BootstrapRepKrzCor(ind.data, iterations = 1000, correlation = F, num.cores = 1)
}
\arguments{
    \item{ind.data}{
        Matrix of residuals or indiviual measurments
    }
    \item{iterations}{
        Number of ressamples to take
    }
    \item{ComparisonFunc}{
        Comparison function for calculated statistic
    }
    \item{StatFunc}{
        Statistic to be calculated
    }
    \item{num.cores}{
        Number of threads to use in computation. 
        Requires doMC library.
    }
    \item{correlation}{
        When using BootstrapRepKrzCor, statistic can be correlation or covariance. If TRUE, uses correlation.
    }

}
\details{
    Bootstraping may be misleading with very small sample sizes. Use
with caution.
}
\value{
    returns the mean repeatability, or mean value of comparisons from
samples to original statistic.
}
\references{
    %% ~put references to the literature/web site here ~
}
\author{
    %%  ~~who you are~~
}
\note{
Diogo Melo, Guilherme Garcia
}
\seealso{
    \code{\link{MonteCarloRep}}, \code{\link{AlphaRep}}
}
\examples{
BootstrapRep(iris[,1:4], iterations = 5,
             ComparisonFunc = function(x, y) RandomSkewers(x, y)[1],
             StatFunc = cov,
             num.cores = 1)

BootstrapRepRandomSkewers(iris[,1:4], 50)

BootstrapRepMantelCor(iris[,1:4], 50)

BootstrapRepKrzCor(iris[,1:4], 50)

BootstrapRepKrzCor(iris[,1:4], 50, TRUE)

#Creating repeatability vector for a list of matrices
mat.list <- RandomMatrix(10, 3, 1, 10)
unlist(lapply(mat.list, BootstrapRepKrzCor), 30, TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bootstrap}
\keyword{repetabilities}% __ONLY ONE__ keyword per line
