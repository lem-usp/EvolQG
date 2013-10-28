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
#'@importFrom phylobase getNode edges reorder
#'@import plyr
#'@examples
#'library(phylobase)
#'data(bird.orders)
#'tree <- as(bird.orders, "phylo4")
#'tip.labels <- tipLabels(tree)
#'mat.list <- RandomMatrix(5, length(tip.labels))
#'names(mat.list) <- tip.labels
#'sample.sizes <- runif(length(tip.labels), 15, 20)
#'AncestralStates(tree, mat.list, sample.sizes)

AncestralStates <- function(tree, tip.data, tip.sample.size = NULL){
  tree <- as(tree, "phylo4")
  tip.names <- names(getNode(tree, type="tip"))
  if(!is.null(tip.sample.size)){
    for(i in 1:length(tip.data)){
      tip.data[[i]] <- tip.data[[i]]*tip.sample.size[i]
    }
  }
  else{
    tip.sample.size <- rep(1, length(tip.data))
  }
  names(tip.sample.size) <- names(tip.data)
  ancestral.stats <- list()
  if(!sum(tip.names %in% names(tip.data)) == length(tip.names)) stop("All tip labels must be in stat list.")
  node.order <- getNode(tree, edges(reorder(tree, order="postorder"))[,2])
  for (node in node.order){
        node.name <- names(getNode(tree, node))
      if(is.na(node.name))
        node.name <- as.character(node)
      descendants.list <- descendants(tree, node)
      ancestral.stats[[node.name]] <- Reduce("+", tip.data[descendants.list])/sum(tip.sample.size[descendants.list])
  }
  return(ancestral.stats)
}
