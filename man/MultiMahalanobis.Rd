\name{MultiMahalanobis}
\alias{MultiMahalanobis}
\title{
    Calculate Mahalonabis distance for many vectors
}
\description{
    Calculates the Mahalanobis distance between a list of species mean,
    using a global covariance matrix
}
\usage{
MultiMahalanobis(mean.list, cov.matrix, num.cores = 1)
}
\arguments{
    \item{mean.list}{
        list of species means being compared
    }
    \item{cov.matrix}{
        covariance matrix definig the metric tensor to be used
    }
    \item{num.cores}{
        Number of threads to use in computation. 
        Requires doMC library.
    }
}

\details{
}
\value{
    returns a distance matrix of species-species comparisons.
}
\references{
    %% ~put references to the literature/web site here ~
}
\author{
    Diogo Melo
}
\note{
    %%  ~~further notes~~
}
\seealso{
    \code{\link{mahalanobis}}, ~~~
}
\examples{
mean.1 <- colMeans(matrix(rnorm(30*10), 30, 10))
mean.2 <- colMeans(matrix(rnorm(30*10), 30, 10))
mean.3 <- colMeans(matrix(rnorm(30*10), 30, 10))
mean.list <- list(mean.1, mean.2, mean.3)

# If cov.matrix is identity, calculated distance is euclidian
euclidian <- MultiMahalanobis(mean.list, diag(rep(1, 10)))
# else, it is not
half.euclidian <- MultiMahalanobis(mean.list, diag(rep(0.5, 10)))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
