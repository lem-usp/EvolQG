\name{TestModularity}
\alias{TestModularity}
\alias{CreateHipotMatrix}
\title{
    Test modularity hipotesis
}
\description{
    Tests modularity hipotesis using cor.matrix matrix and trait groupings
}
\usage{
    TestModularity(cor.matrix, modularity.hipot, iterations = 100)
}
\arguments{
    \item{cor.matrix}{
        Correlation matrix
    }
    \item{modularity.hipot}{
        modularity.hipot: Matrix of hipotesis.
        Each line represents a trait and each column a module.
        if modularity.hipot[i,j] == 1, trait i is in module j.
    }
    \item{iterations}{
        Number of iterations, to be passed to MantelCor
    }
}
\details{
    %%  ~~ If necessary, more details than the description above ~~
}
\value{
    Returns mantel correlation and associated probability for each modularity hipotesis, along with AVG+, AVG-, AVG Ratio for each module.
    A total hipotesis combining all hipotesis is also tested.
}
\references{
    %% ~put references to the literature/web site here ~
}
\author{
    Diogo Melo, Guilherme Garcia
}
\note{
    %%  ~~further notes~~
}

\seealso{
    \code{\link{MantelCor}}
}
\examples{
    rand.hipots <- matrix(sample(c(1, 0), 30, replace=T), 10, 3)
    hip.array <- CreateHipotMatrix(rand.hipots)
    expect_that(sum(laply(hip.array, isSymmetric)), equals(4))
    expect_that(hip.array, is_a("list"))
    expect_that(dim(hip.array[[1]]), equals(c(10, 10)))
    mod.test <- TestModularity(cor.matrix, rand.hipots)
}
\keyword{mantel}
\keyword{modularity}
