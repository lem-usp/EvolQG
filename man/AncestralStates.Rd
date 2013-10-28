\name{AncestralStates}
\alias{AncestralStates}
\title{Calculates ancestral states of some statistic}
\usage{
  AncestralStates(tree, tip.data, tip.sample.size = NULL)
}
\arguments{
  \item{tree}{phylogenetic tree}

  \item{tip.data}{list of tip nodes data}

  \item{tip.sample.size}{vector of tip nodes sample sizes}
}
\value{
  list with calculated ancestral states, using labels or
  numbers from tree
}
\description{
  Calculates weighted average of some statistic along a
  phylogeny
}
\examples{
library(phylobase)
data(bird.orders)
tree <- as(bird.orders, "phylo4")
tip.labels <- tipLabels(tree)
mat.list <- RandomMatrix(5, length(tip.labels))
names(mat.list) <- tip.labels
sample.sizes <- runif(length(tip.labels), 15, 20)
AncestralStates(tree, mat.list, sample.sizes)
}

