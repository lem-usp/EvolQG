\name{CalculateMatrix}
\alias{CalculateMatrix}
\title{Calculate Covariance Matrix from LinearModel}
\usage{
  CalculateMatrix(linear.m)
}
\arguments{
  \item{linear.m}{Linear model adjusted for original data.}
}
\description{
  Calculates covariance matrix using the maximum likelihood
  estimator and the model residuals.
}
\examples{
data(iris)
options(contrasts=c("contr.sum","contr.poly"))
iris.lm = lm(as.matrix(iris[,1:4])~iris[,5])
cov.matrix <- CalculateMatrix(iris.lm)

#To obtain a corrlation matrix, use:
cor.matrix <- cov2cor(cov.matrix)
}
\author{
  Diogo Melo, Fabio Machado
}
\references{
  https://github.com/lem-usp/Morphometrics/wiki/
}

