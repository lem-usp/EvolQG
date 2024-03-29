% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MINT.R
\name{MINT}
\alias{MINT}
\alias{JackKnifeMINT}
\title{Modularity and integration analysis tool}
\usage{
MINT(
  c.matrix,
  modularity.hypot,
  significance = FALSE,
  sample.size = NULL,
  iterations = 1000
)

JackKnifeMINT(
  ind.data,
  modularity.hypot,
  n = 1000,
  leave.out = floor(dim(ind.data)[1]/10),
  ...
)
}
\arguments{
\item{c.matrix}{Correlation or covariance matrix}

\item{modularity.hypot}{Matrix of hypothesis. Each line represents a trait and each column a module.
if modularity.hypot[i,j] == 1, trait i is in module j.}

\item{significance}{Logical. Indicates if goodness of fit test should be performed.}

\item{sample.size}{sample size in goodness of fit simulations via MonteCarlo}

\item{iterations}{number os goodness of fit simulations}

\item{ind.data}{Matrix of residuals or individual measurements}

\item{n}{number of jackknife samples}

\item{leave.out}{number of individuals to be left out of each jackknife, default is 10\%}

\item{...}{additional arguments to be passed to raply for the jackknife}
}
\value{
Dataframe with ranked hypothesis, ordered by the corrected gamma value
Jackknife will return the best hypothesis for each sample.
}
\description{
Combines and compares many modularity hypothesis to a covariance matrix. Comparison values are
adjusted to the number of zeros in the hypothesis using a linear regression. Best hypothesis can
be assessed using a jack-knife procedure.
}
\note{
Hypothesis can be named as column names, and these will be used to make labels in the 
output.
}
\examples{
# Creating a modular matrix:
modules = matrix(c(rep(c(1, 0, 0), each = 5),
                 rep(c(0, 1, 0), each = 5),
                 rep(c(0, 0, 1), each = 5)), 15)
                 
cor.hypot = CreateHypotMatrix(modules)[[4]]
hypot.mask = matrix(as.logical(cor.hypot), 15, 15)
mod.cor = matrix(NA, 15, 15)
mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.1, 0.2) # between-modules
diag(mod.cor) = 1
mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric

# True hypothesis and a bunch of random ones.
hypothetical.modules = cbind(modules, matrix(sample(c(1, 0), 4*15, replace=TRUE), 15, 4))

# if hypothesis columns are not named they are assigned numbers
colnames(hypothetical.modules) <- letters[1:7]
 
MINT(mod.cor, hypothetical.modules)

random_var = runif(15, 1, 10)
mod.data = mvtnorm::rmvnorm(100, sigma = sqrt(outer(random_var, random_var)) * mod.cor)
out_jack = JackKnifeMINT(mod.data, hypothetical.modules, n = 50)

library(ggplot2)

ggplot(out_jack, aes(rank, corrected.gamma)) + geom_point() + 
       geom_errorbar(aes(ymin = lower.corrected, ymax = upper.corrected))
}
\references{
Marquez, E.J. 2008. A statistical framework for testing modularity in multidimensional data. 
Evolution 62:2688-2708. 

Parsons, K.J., Marquez, E.J., Albertson, R.C. 2012. Constraint and opportunity: the genetic 
basis and evolution of modularity in the cichlid mandible. The American Naturalist 179:64-78.

http://www-personal.umich.edu/~emarquez/morph/doc/mint_man.pdf
}
