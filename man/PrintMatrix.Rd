\name{PrintMatrix}
\alias{PrintMatrix}
\alias{PrintMatrix.default}
\alias{PrintMatrix.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print Matrix to file
}
\description{
Print a matrix or a list of matrices to file
}
\usage{
PrintMatrix(x, ...)
\method{PrintMatrix}{default}(x, out.file = './matrix.csv', ...)
\method{PrintMatrix}{list}(x, out.file = './matrix.csv', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix or list of matrices
}
  \item{out.file}{
Output file
}
  \item{...}{
Aditional parameters
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    Prints coma separated matrices, with labels
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Diogo Melo
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
m.list <- RandomMatrix(10, 4)
PrintMatrix(m.list)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
