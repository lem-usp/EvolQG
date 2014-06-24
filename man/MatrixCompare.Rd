\name{MatrixCompare}
\alias{MatrixCompare}
\title{Matrix Compare}
\usage{
MatrixCompare(cov.x, cov.y)
}
\arguments{
  \item{cov.x}{covariance or correlation matrix}

  \item{cov.y}{covariance or correlation matrix}
}
\value{
data.frame of comparisons
}
\description{
Compare two matrices using all available methods. Currently
RandomSkewers, MantelCor and KrzCor
}
\examples{
cov.x = RandomMatrix(10, 1, 1, 10)
cov.y = RandomMatrix(10, 1, 10, 20)
MatrixCompare(cov.x, cov.y)
}

