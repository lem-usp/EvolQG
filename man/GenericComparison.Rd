\name{Comparison Map}
\alias{ComparisonMap}
\alias{SingleComparisonMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Generic Comparison Map functions for creating parallel list methods
}
\description{
    %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ComparisonMap(matrix.list, MatrixCompFunc, repeat.vector = NULL,
              num.cores = 1)
SingleComparisonMap(matrix.list, y.mat, MatrixCompFunc, num.cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{matrix.list}{
        list of matrices being compared
    }
    \item{y.mat}{
        single matrix being compared to list
    }
    \item{MatrixCompFunc}{
        Function used to compare pair of matrices
    }
    \item{repeat.vector}{
        Vector of repeatabilities for correlation correction.
    }
    \item{num.cores}{
        Number of threads to use in computation. 
        Requires doMC library.
    }
}
\details{
    %%  ~~ If necessary, more details than the description above ~~
}
\value{
    Matrix of comparisons, matrix of probabilities.

    data.frame of comparisons and probabilities.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
