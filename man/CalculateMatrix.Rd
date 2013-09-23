\name{CalculateMatrix}
\alias{CalculateMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Covariance Matrix from LinearModel
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CalculateMatrix(linear.m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linear.m}{
%%     ~~Describe \code{linear.m} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (linear.m) 
{
    cov.matrix = var(linear.m$residuals) * ((dim(linear.m$residuals)[1] - 
        1)/linear.m$df.residual)
    return(cov.matrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
