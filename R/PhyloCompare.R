#'Compares sister groups
#'
#'Calculates the comparison of some statistic between sister groups along a phylogeny
#'
#'@param tree phylogenetic tree
#'@param node.data list of node data
#'@param CompareFunc comparison function
#'@return list with calculated comparisons for each node, using labels or numbers from tree
#'@note Phylogeny must be fully resolved
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
#'phyl.state <- AncestralStates(tree, mat.list, sample.sizes)
#'PhyloCompare(tree, phyl.state)
#'
PhyloCompare <- function(tree, node.data, CompareFunc = RandomSkewers){
  if(is.null(tree$node.label)){
    node.names <- tree$tip.label
  } else{
    node.names <- c(tree$tip.label, tree$node.label)
  }
  # if(!all(tree$tip.label %in% names(node.data))) stop("All tree node labels must be in node.data.")
  node.order <- unique(reorder(tree, "postorder")$edge[,1])
  phylo.comparisons <- list()
  for (node in node.order){
    if(is.na(node.names[node])){
      node.names[node] <- as.character(node)
    }
    current.nodes <- tree$edge[which(tree$edge[,1]==node),2]
    phylo.comparisons[[node.names[node]]] <- CompareFunc(node.data[[current.nodes[1]]], node.data[[current.nodes[2]]])
  }
  phylo.comparisons = ldply(phylo.comparisons, .id = 'node')
  return(phylo.comparisons)
}
