\name{Trait Repeatability}
\alias{CalcRepeatability}
\title{
    Parametric per trait Repeatabilities
}
\description{
    Estimates the variance in the sample not due to measurment error
}
\usage{
CalcRepeatability(ID, ind.data)
}
\arguments{
    \item{ID}{
        Individual ID labels.
    }
    \item{ind.data}{
        matrix of residuals or individual measurments
    }
}
\details{
    Requires at least two observations per individual
}
\value{
    Vector of Repeatabilities per trait
}
\references{
    Lessels, C. M., & Boag, P. T. (1987).
    Unrepeatable repeatabilities: a common mistake.
    The Auk, 2(January), 116¿121.
}
\author{
    Guilherme Garcia
}
\note{
    %%  ~~further notes~~
}

\seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
num.ind = length(iris[,1])
ID = rep(1:num.ind, 2)
ind.data = rbind(iris[,1:4], iris[,1:4]+array(rnorm(num.ind*4, 0, 0.1), dim(iris[,1:4])))
CalcRepeatability(ID, ind.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Repeatability}
