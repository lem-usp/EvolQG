% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rarefaction.R
\name{Rarefaction}
\alias{Rarefaction}
\title{Rarefaction analysis via resampling}
\usage{
Rarefaction(
  ind.data,
  ComparisonFunc,
  ...,
  num.reps = 10,
  correlation = FALSE,
  replace = FALSE,
  parallel = FALSE
)
}
\arguments{
\item{ind.data}{Matrix of residuals or individual measurments}

\item{ComparisonFunc}{comparison function}

\item{...}{Additional arguments passed to ComparisonFunc}

\item{num.reps}{number of populations sampled per sample size}

\item{correlation}{If TRUE, correlation matrix is used, else covariance matrix. MantelCor always uses correlation matrix.}

\item{replace}{If true, samples are taken with replacement}

\item{parallel}{if TRUE computations are done in parallel. Some foreach back-end must be registered, like doParallel or doMC.}
}
\value{
returns the mean value of comparisons from samples to original statistic, for all sample sizes.
}
\description{
Calculates the repeatability of a statistic of the data, such as
correlation or covariance matrix, via bootstrap resampling with
varying sample sizes, from 2 to the size of the original data.
}
\details{
Samples of various sizes, with replacement, are taken from the full population, a statistic calculated
and compared to the full population statistic.

A specialized plotting function displays the results in publication quality.

Bootstraping may be misleading with very small sample sizes. Use with caution if original sample sizes are small.
}
\examples{
ind.data <- iris[1:50,1:4]

results.RS <- Rarefaction(ind.data, RandomSkewers, num.reps = 5)
#' #Easy parsing of results
library(reshape2)
melt(results.RS)

# or :
\donttest{
results.Mantel <- Rarefaction(ind.data, MatrixCor, correlation = TRUE, num.reps = 5)
results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5)
results.PCA <- Rarefaction(ind.data, PCAsimilarity, num.reps = 5)
}

\dontrun{
#Multiple threads can be used with some foreach backend library, like doMC or doParallel
library(doMC)
registerDoMC(cores = 2)
results.KrzCov <- Rarefaction(ind.data, KrzCor, num.reps = 5, parallel = TRUE)
}
}
\seealso{
\code{\link{BootstrapRep}}
}
\author{
Diogo Melo, Guilherme Garcia
}
\keyword{bootstap}
\keyword{rarefaction}
\keyword{repeatability}
