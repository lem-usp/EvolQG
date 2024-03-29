% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/localShapeVariables.R
\name{LocalShapeVariables}
\alias{LocalShapeVariables}
\title{Local Shape Variables}
\usage{
LocalShapeVariables(
  gpa = NULL,
  cs = NULL,
  landmarks = NULL,
  tesselation,
  run_parallel = FALSE
)
}
\arguments{
\item{gpa}{Procustes aligned landmarks.}

\item{cs}{Centoid sizes}

\item{landmarks}{unaligned landmarks. Ignored if both gpa and cs are passed.}

\item{tesselation}{matrix of rows of the landmarks. The centroid of each row
is used to mark the position of the jacobians}

\item{run_parallel}{Logical. If computation should be paralleled. Use with
caution, can make things worse. Requires that at parallel back-end like doMC
be registered}
}
\value{
List with TPS functions, jacobian matrices, local shape variables, mean shape, centroid sizes and individual IDs
}
\description{
Calculates the local shape variables of a set of landmarks using the sequence:
- TPS transform between all shapes and the mean shape
- Jacobian of the TPS transforms at the centroid of rows of the landmarks in the tesselation argument
- Mean center the Jacobians using the Karcher Mean
- Take the determinant of the centered jacobians
}
\author{
Guilherme Garcia

Diogo Melo
}
