\name{CalculateMatrix}
\alias{CalculateMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Calculate Covariance Matrix from LinearModel
}
\description{
    Calculates covariance matrix using the maximum likelihood estimator and the model residuals.
}
\usage{
CalculateMatrix(linear.m)
}
\arguments{
    \item{linear.m}{
        Linear model adjusted for original data.
    }
}
\details{
    %%  ~~ If necessary, more details than the description above ~~
}
\value{
    Estimated covariance matrix.
}
\references{
    https://github.com/lem-usp/lem-usp-R/wiki/
}
\author{
    Diogo Melo, Fabio Machado
}
\note{
    %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(iris)
options(contrasts=c("contr.sum","contr.poly"))
iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
cov.matrix <- CalculateMatrix(iris.lm)

#To obtain a corrlation matrix, use:
cor.matrix <- cov2cor(cov.matrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{covariancematrix}
