\name{PlotSRD}
\alias{PlotSRD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotSRD(output, matrix.label = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
%%     ~~Describe \code{output} here~~
}
  \item{matrix.label}{
%%     ~~Describe \code{matrix.label} here~~
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
function (output, matrix.label = "") 
{
    layout(array(c(1, 1, 2, 2), c(2, 2)))
    par(mar = c(4, 4, 8, 0.4))
    mean.r2 <- output$out[, 1]
    low.r2 <- output$out[, 2]
    up.r2 <- output$out[, 3]
    c.mean.r2 <- output$out[, 5]
    c.sd.r2 <- output$out[, 6]
    if (is.null(rownames(output$out))) 
        dists <- 1:length(mean.r2)
    else dists <- rownames(output$out)
    pc.pch <- output$model$code + 20
    plot(mean.r2, type = "p", lty = 2, pch = pc.pch, ylab = "", 
        xlab = "", xaxt = "n", ylim = c(-1, 1))
    for (i in 1:length(mean.r2)) {
        abline(v = i, lty = 1, col = rgb(0.8, 0.8, 0.8))
    }
    arrows(x0 = 1:length(mean.r2), y0 = low.r2, y1 = up.r2, angle = 90, 
        length = 0.05, code = 3)
    abline(h = mean(mean.r2), lty = 3)
    axis(3, 1:length(mean.r2), dists, las = 2, cex.axis = 1.3)
    mtext(side = 2, at = mean(mean.r2), text = round(mean(mean.r2), 
        2), las = 2, cex = 1.8)
    par(mar = c(4, 0, 8, 4.6))
    output$model$code[output$model$code == 1] <- 1.2
    output$model$code[output$model$code == 0] <- 1
    pc.cex <- output$model$code
    plot(c.sd.r2 ~ c.mean.r2, pch = pc.pch, xlab = "", ylab = "", 
        yaxt = "n", main = matrix.label, cex.main = 3, cex.axis = 1.3)
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
    text(c.mean.r2, c.sd.r2, labels = dists, pos = 4, cex = pc.cex)
    axis(4, las = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
