\name{plot.group.evol}
\alias{plot.group.evol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.group.evol(evol.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{evol.list}{
%%     ~~Describe \code{evol.list} here~~
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
function (evol.list) 
{
    upper.crit <- function(vec, val = 0.95) {
        n.obs <- length(vec)
        vec <- sort(vec)
        crit.pos <- round(val * n.obs)
        return(vec[crit.pos])
    }
    plot.single <- function(evol.list, which) {
        n.obj <- length(evol.list)
        extract.dist <- function(element, which) return(element$dist[, 
            which]/element$max.val[which])
        to.plot <- sapply(evol.list, extract.dist, which)
        colnames(to.plot) <- names(evol.list)
        boxplot(to.plot, cex = 0.2, boxwex = 0.5, border = rgb(0.6, 
            0.6, 0.6))
        extract.mod <- function(element, which) return(element$mod[, 
            which]/element$max.val[which])
        extract.mod.names <- function(element) return(rownames(element$mod))
        mod.values <- lapply(evol.list, extract.mod, which)
        mod.names <- lapply(evol.list, extract.mod.names)
        crits <- apply(to.plot, 2, upper.crit)
        for (i in 1:n.obj) {
            segments(x0 = i - 0.25, x1 = i + 0.25, y0 = crits[i], 
                col = "red", lwd = 2)
            segments(x0 = i - 0.25, x1 = i + 0.25, y0 = mod.values[[i]], 
                col = "black", lty = i)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
