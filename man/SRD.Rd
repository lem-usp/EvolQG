\name{SRD}
\alias{PlotSRD}
\alias{SRD}
\alias{SRD.default}
\alias{SRD.list}
\title{Compare matrices via Selection Response Decomposition}
\usage{
  SRD(cov.x, cov.y, ...)

  \method{SRD}{default} (cov.x, cov.y, iterations = 1000,
    ...)

  \method{SRD}{list} (cov.x, cov.y = NULL,
    iterations = 1000, num.cores = 1, ...)

  PlotSRD(output, matrix.label = "")
}
\arguments{
  \item{cov.x}{Covariance matrix being compared. cov.x can
  be a matrix or a list.}

  \item{cov.y}{Covariance matrix being compared. Ignored if
  cov.x is a list.}

  \item{iterations}{Number of random vectors used in
  comparison}

  \item{num.cores}{If list is passed, number of threads to
  use in computation. Requires doMC library.}

  \item{...}{aditional parameters passed to other methods}

  \item{output}{Output form SRD function, used in ploting}

  \item{matrix.label}{Plot label}
}
\value{
  out List of SRD scores mean, confidence intervals,
  standard deviation, centerd mean e centerd standard
  deviation

  pc1 scored along the pc1 of the mean/SD correlation
  matrix

  model List of linear model results from mean/SD
  correlation. Quantiles, interval and divergent traits

  cormat Mean/Sd model correlation matrix
}
\description{
  Based on Random Skewers tecnique, selection response
  vectors are expanded in direct and indirect componenet by
  trait and compared via vector correlations.
}
\details{
  Output can be ploted using PlotSRD function
}
\note{
  If input is a list, output is a symmetric list array with
  pairwise comparisons.
}
\examples{
cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
colnames(cov.matrix.1) <- colnames(cov.matrix.2) <- sample(letters, 10)
rownames(cov.matrix.1) <- rownames(cov.matrix.2) <- colnames(cov.matrix.1)
srd.out <- SRD(cov.matrix.1, cov.matrix.2)

#lists
m.list <- RandomMatrix(10, 4)
SRD(m.list)

#divergent traits
colnames(cov.matrix.1)[as.logical(srd.out$model$code)]

#Plot
PlotSRD(srd.out)
}
\author{
  Diogo Melo, Guilherme Garcia
}
\references{
  Marroig, G., Melo, D., Porto, A., Sebastiao, H., and
  Garcia, G. (2011). Selection Response Decomposition
  (SRD): A New Tool for Dissecting Differences and
  Similarities Between Matrices. Evolutionary Biology,
  38(2), 225-241. doi:10.1007/s11692-010-9107-2
}
\seealso{
  \code{\link{RandomSkewers}}
}

