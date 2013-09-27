\name{SRD}
\alias{SRD}
\alias{PlotSRD}
\title{
    Compare matrices via Selection Response Decomposition
}
\description{
    Based on Random Skewers tecnique, selection response vectors are
    expanded in direct and indirect componenet by trait and compared via
    vector correlations.
}
\usage{
SRD(cov.matrix.1, cov.matrix.2, iterations = 1000)
PlotSRD (output, matrix.label = "")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cov.matrix.1, cov.matrix.2}{
        Covariance matrix being compared
    }
    \item{iterations}{
        Number of random vectors used in comparison
    }
    \item{output}{
        Output form SRD function, used in ploting
    }
    \item{matrix.label}{
        Plot label
    }
}
\details{
    Output can be ploted using PlotSRD function
}
\value{
    \item{out}{List of SRD scores mean, confidence intervals, standard
        deviation, centerd mean e centerd standard deviation}
    \item{pc1}{scored along the pc1 of the mean/SD correlation matrix}
    \item{model}{List of linear model results from mean/SD correlation. Quantiles, interval and divergent traits}
    \item{cormat}{Mean/Sd model correlation matrix}

}
\references{
    Marroig, G., Melo, D., Porto, A., Sebastiao, H., and Garcia, G.
    (2011). Selection Response Decomposition (SRD): A New Tool for
    Dissecting Differences and Similarities Between Matrices. Evolutionary
    Biology, 38(2), 225-241. doi:10.1007/s11692-010-9107-2
}
\author{
    Diogo Melo, Guilherme Garcia
}
\note{
    %%  ~~further notes~~
}
\seealso{
    \code{\link{RandomSkewers}}
}
\examples{
cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
colnames(cov.matrix.1) <- colnames(cov.matrix.2) <- sample(letters, 10)
rownames(cov.matrix.1) <- rownames(cov.matrix.2) <- colnames(cov.matrix.1)
srd.out <- SRD(cov.matrix.1, cov.matrix.2)


#divergent traits
colnames(cov.matrix.1)[as.logical(srd.out$model$code)]

#Plot
PlotSRD(srd.out)

}
\keyword{SRD}
\keyword{RandomSkewers}
\keyword{selectionresponsedecomposition}% __ONLY ONE__ keyword per line
