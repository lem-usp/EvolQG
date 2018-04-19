#' Mantel test with phylogenetic permutations
#' 
#' Performs a matrix correlation with significance given by a phylogenetic Mantel Test. 
#' Pairs of rowns and columns are permuted with probability proportional to their phylogenetic distance.
#' 
#' @param tree phylogenetic tree. Tip labels must match names in input matrices
#' @param matrix.1 pair-wise comparison/distance matrix
#' @param matrix.2 pair-wise comparison/distance matrix
#' @param ... aditional parameters, currently ignored
#' @param permutations Number of permutations used in significance calculation
#' @param ComparisonFunc comparison function, default is MatrixCor
#' @param k determines the influence of the phylogeny. 1 is strong influence, and larger values converge to a traditional mantel test.
#' @return returns a vector with the comparison value and the proportion of times the observed comparison is smaller than the correlations from the permutations. 
#' @export
#' @note This method should only be used when there is no option other than representing data as pair-wise. It suffers from low power, and alternatives should be used when available.
#' @references Harmon, L. J., & Glor, R. E. (2010). Poor statistical performance of the Mantel test in phylogenetic comparative analyses. Evolution, 64(7), 2173-2178.
#' @references Lapointe, F. J., & Garland, Jr, T. (2001). A generalized permutation model for the analysis of cross-species data. Journal of Classification, 18(1), 109-127.
#' @author Diogo Melo, adapted from Harmon & Glor 2010
#' @examples 
#' data(dentus)
#' data(dentus.tree)
#' tree = dentus.tree
#' cor.matrices = dlply(dentus, .(species), function(x) cor(x[-5]))
#' comparisons = MatrixCor(cor.matrices)
#' 
#' sp.means = dlply(dentus, .(species), function(x) colMeans(x[-5]))
#' mh.dist = MultiMahalanobis(means = sp.means, cov.matrix = PhyloW(dentus.tree, cor.matrices)$'6')
#' PhyloMantel(dentus.tree, comparisons, mh.dist, k = 10000)
#' 
#' #similar to MantelCor for large k:
#' \dontrun{
#' PhyloMantel(dentus.tree, comparisons, mh.dist, k = 10000)
#' MantelCor(comparisons, mh.dist)
#' }
PhyloMantel <- function(tree, matrix.1, matrix.2, ..., 
                        permutations=1000, 
                        ComparisonFunc = function(x, y) cor(x[lower.tri(x)], y[lower.tri(y)]), 
                        k = 1){ 
  corr <- ComparisonFunc(matrix.1, matrix.2) 
  nullstats <- raply(permutations, 
                     cor(matrix.1[lower.tri(matrix.1)], 
                         permPhylo(tree, matrix.2, k = k)[lower.tri(matrix.1)])) 
  prob <- sum(nullstats >= corr)/(permutations)
  c("Rsquared" = corr, "Probability" = prob)
}

# Calculates probabilities for phylogenetic permutations; from Lapointe and Garland 2001
#' @importFrom ape cophenetic.phylo
phyloProb<-function(tree, k) {
  phylo_dist <- cophenetic.phylo(tree)
  scaled_phylo_dist <- phylo_dist/max(phylo_dist)
  s <- k-scaled_phylo_dist
  s/rowSums(s)
}

# Permutes species according to phylogentic tree; returns tip names in permuted order
phyloPermute<-function(tree, k) {
  p<-phyloProb(tree, k)
  tt<-rownames(p)
  nsp<-length(tt)
  order<-sample(1:nsp, replace=F)
  ttnew<-character(nsp)
  cpm<-p
  for(j in order[-nsp]) {
    cpm<-cpm/rowSums(cpm)
    rr<-which(rownames(cpm)==tt[j])
    pp<-cpm[rr,]
    s2<-sample(names(pp), size=1, replace=T, prob=pp)
    slot<-which(tt==s2)
    rc<-which(colnames(cpm)==s2)
    ttnew[slot]<-tt[j]
    cpm<-cpm[-rr,-rc]
  }
  ttnew[which(ttnew=="")]<-tt[order[nsp]]
  
  ttnew
}	

permPhylo <- function(tree, matrix.1, k) 
{ 
  s <- phyloPermute(tree, k = k) 
  trans <- match(s, rownames(matrix.1))
  matrix.1[trans,trans] 
}