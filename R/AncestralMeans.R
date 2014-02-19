#'Calculates ancestral states of some statistic
#'
#'Calculates weighted average of some statistic along a phylogeny
#'
#'@param tree phylogenetic tree
#'@param tip.data list of tip nodes data
#'@param tip.sample.size vector of tip nodes sample sizes
#'
#'@return list with calculated ancestral states, using labels or numbers from tree
#'@export
#'@importFrom ape reorder.phylo
#'@import plyr
#'@examples
#'library(ape)
#'data(bird.orders)
#'tree <- bird.orders
#'mat.list <- RandomMatrix(5, length(tree$tip.label))
#'names(mat.list) <- tree$tip.label
#'sample.sizes <- runif(length(tree$tip.label), 15, 20)
#'AncestralStates(tree, mat.list, sample.sizes)

AncestralStates <- function(tree, tip.data, tip.sample.size = NULL){
  num.nodes = length(tree$tip.label)
  if(is.null(tree$node.label))
    node.names <- tree$tip.label 
  else
    node.names <- c(tree$tip.label, tree$node.label)
  if(is.null(tip.sample.size))
    tip.sample.size <- rep(1, length(tip.data))
  names(tip.sample.size) <- names(tip.data)
  tip.sample.size <- as.list(tip.sample.size)
  ancestral.stats <- tip.data
  if(!all(tree$tip.label %in% names(tip.data))) stop("All tip labels must be in stat list.")
  node.order <- unique(reorder(tree, "postorder")$edge[,1])
  for (node in node.order){
        node.name <- node.names[node]
      if(is.na(node.name)){
        node.names[node] <- as.character(node)
      }
      descendants.list <- node.names[tree$edge[which(tree$edge[,1]==node),2]]
      ponderados <- Map('*', ancestral.stats[descendants.list], tip.sample.size[descendants.list])
      tip.sample.size[[node.names[node]]] <- Reduce("+", tip.sample.size[descendants.list])
      ancestral.stats[[node.names[node]]] <- Reduce("+", ponderados)/
                                              tip.sample.size[[node.names[node]]]
  }
  return(ancestral.stats)
}