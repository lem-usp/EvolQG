\name{CalcR2CvCorrected}
\alias{CalcR2CvCorrected}
\alias{CalcR2CvCorrected.default}
\alias{CalcR2CvCorrected.lm}
\title{
    Corrected integration value
}
\description{
    Calculates the Young correction for integration, using bootstrap
}
\usage{
CalcR2CvCorrected(ind.data, ...)
\method{CalcR2CvCorrected}{default}(ind.data, iterations = 1000, num.cores = 1, ...)
\method{CalcR2CvCorrected}{lm}(ind.data, iterations = 1000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{ind.data}{
        Matrix of residuals or indiviual measurments, or ajusted linear model
    }
    \item{iterations}{
        Number of ressamples to take
    }
    \item{num.cores}{
        Number of threads to use in computation.
        Requires doMC library.
    }
    \item{...}{aditional arguments passed to other methods}
}
\details{
    %%  ~~ If necessary, more details than the description above ~~
}
\value{
    data.frame with mean squared correlation, mean eigen value coeficient of variation and mean coeficient of variation for all sampled populations
}
\references{
    Young, N. M., Wagner, G. P., and Hallgrimsson, B. (2010).
Development and the evolvability of human limbs. Proceedings of the
National Academy of Sciences of the United States of America, 107(8),
3400-5. doi:10.1073/pnas.0911856107
}
\author{
    Diogo Melo, Guilherme Garcia
}
\note{
    %%  ~~further notes~~
}

\seealso{
    \code{\link{MeanMatrixStatistics}}, \code{\link{CalcR2}}
}
\examples{
r2.dist = CalcR2CvCorrected(iris[,1:4])
plot(r2.dist[,1], r2.dist[,3])
plot(r2.dist[,2], r2.dist[,3])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{correlations}
\keyword{integrations}
