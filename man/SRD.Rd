\name{SRD}
\alias{SRD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SRD(cov.matrix.1, cov.matrix.2, nsk = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.matrix.1}{
%%     ~~Describe \code{cov.matrix.1} here~~
}
  \item{cov.matrix.2}{
%%     ~~Describe \code{cov.matrix.2} here~~
}
  \item{nsk}{
%%     ~~Describe \code{nsk} here~~
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
function (cov.matrix.1, cov.matrix.2, nsk = 1000) 
{
    size <- dim(cov.matrix.1)[1]
    r2s <- array(0, c(size, nsk))
    beta <- apply(array(rnorm(size * nsk, mean = 0, sd = 1), 
        c(size, nsk)), 2, Normalize)
    for (I in 1:nsk) {
        beta.matrix <- diag(beta[, I])
        dz1 <- apply(cov.matrix.1 \%*\% beta.matrix, 1, Normalize)
        dz2 <- apply(cov.matrix.2 \%*\% beta.matrix, 1, Normalize)
        r2s[, I] <- colSums(dz1 * dz2)
    }
    mean.r2 <- apply(r2s, 1, mean)
    sample.conf <- function(x, lower = TRUE) {
        ox <- x[order(x)]
        lox <- length(ox)
        if (lower) 
            crit <- round(0.025 * lox)
        else crit <- round(0.975 * lox)
        return(ox[crit])
    }
    low.r2 <- apply(r2s, 1, sample.conf, lower = TRUE)
    up.r2 <- apply(r2s, 1, sample.conf, lower = FALSE)
    sd.r2 <- apply(r2s, 1, sd)
    cmean.r2 <- scale(mean.r2, scale = FALSE)
    csd.r2 <- scale(sd.r2, scale = FALSE)
    cent <- cbind(cmean.r2, csd.r2)
    pca.cent <- princomp(cent, cor = TRUE, scores = TRUE)
    pc1 <- pca.cent$scores[, 1]
    if (pca.cent$loadings[1, 1] < 0) 
        pc1 <- -pc1
    pc1 <- pc1/sd(pc1)
    pc1.quant <- quantile(pc1, probs = c(1, 5, 10, 20, 25, 30, 
        40, 50, 60, 70, 75, 80, 90, 95, 99)/100, names = FALSE)
    pc1.int <- -1.96/sqrt(length(pc1))
    pc1.sig <- ifelse(pc1 < pc1.int, 1, 0)
    model <- list(quantiles = pc1.quant, interval = pc1.int, 
        code = pc1.sig)
    output <- cbind(mean.r2, low.r2, up.r2, sd.r2, cmean.r2, 
        csd.r2)
    colnames(output) <- c("ARC", "IC-", "IC+", "SD", "CMEAN", 
        "CSD")
    rownames(output) <- rownames(cov.matrix.1)
    return(list(out = output, pc1 = pc1, model = model, cormat = cor(t(r2s))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
