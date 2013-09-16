\name{plot.mod.evol}
\alias{plot.mod.evol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.mod.evol(evo.out, new.dev = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{evo.out}{
%%     ~~Describe \code{evo.out} here~~
}
  \item{new.dev}{
%%     ~~Describe \code{new.dev} here~~
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
function (evo.out, new.dev = TRUE) 
{
    require(MASS)
    with(evo.out, {
        n.hip <- nrow(mod)
        n.sk <- nrow(dist)
        if (new.dev) 
            par(mfrow = c(2, 3))
        out <- array(0, dim(mod))
        for (j in 1:6) {
            dist[, j] <- dist[, j]/max.val[j]
            mod[, j] <- mod[, j]/max.val[j]
            x.lim <- c(ifelse(min(mod[, j]) < min(dist[, j]), 
                min(mod[, j]), min(dist[, j])), ifelse(max(mod[, 
                j]) > max(dist[, j]), max(mod[, j]), max(dist[, 
                j])))
            truehist(dist[, j], prob = TRUE, border = "grey", 
                col = "white", xlim = x.lim, xlab = colnames(mod)[j], 
                main = "")
            lines(density(dist[, j]))
            abline(v = sort(dist[, j])[round(0.975 * n.sk)], 
                lty = 2, col = "red")
            abline(v = sort(dist[, j])[round(0.025 * n.sk)], 
                lty = 2, col = "red")
            for (i in 1:n.hip) {
                abline(v = mod[i, j], lty = 3, col = rgb(0.2, 
                  0.2, 0.2))
                mtext(side = 3, at = mod[i, j], text = rownames(mod)[i], 
                  las = 3, cex = 0.7)
                out[i, j] <- sum(mod[i, j] < dist[, j])/n.sk
            }
        }
        dimnames(out) <- dimnames(mod)
        return(out)
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
