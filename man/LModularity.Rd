\name{LModularity}
\alias{LModularity}
\title{L Modularity}
\usage{
  LModularity(cor.matrix)
}
\arguments{
  \item{cor.matrix}{correlation matrix}
}
\value{
  List with L-Modulariy value and trait partition
}
\description{
  Calculates the L-Modularity and a partition of traits
}
\examples{
cor.matrix = RandomMatrix(10)
LModularity(cor.matrix)
}

