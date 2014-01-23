\name{PrintMatrix}
\alias{PrintMatrix}
\alias{PrintMatrix.default}
\alias{PrintMatrix.list}
\title{Print Matrix to file}
\usage{
PrintMatrix(x, ...)

\method{PrintMatrix}{default}(x, out.file = "./matrix.csv", ...)

\method{PrintMatrix}{list}(x, out.file = "./matrix.csv", ...)
}
\arguments{
  \item{x}{Matrix or list of matrices}

  \item{out.file}{Output file}

  \item{...}{Aditional parameters}
}
\value{
Prints coma separated matrices, with labels
}
\description{
Print a matrix or a list of matrices to file
}
\examples{
m.list <- RandomMatrix(10, 4)
PrintMatrix(m.list)
}
\author{
Diogo Melo
}

