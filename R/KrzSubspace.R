#' Krzanowski common subspaces analysis
#' 
#' Calculates the subspace most similar across a set of covariance matrices.
#' 
#' @param cov.matrices list of covariance matrices
#' @param k number of dimentions to be retained in calculating the subspace
#' @return H shared space matrix
#' @return k_eVals_H eigen values for shared space matrix, maximum value for each is the number of matrices, representing a fully shared direction
#' @return k_eVecs_H eigen vectors of shared space matrix
#' @return angles between each population subspace and each eigen vector of shared space matrix
#' @note can be used to inplement the Baysean comparsion from Aguirre et al. 2014
#' @references Aguirre, J. D., E. Hine, K. McGuigan, and M. W. Blows. "Comparing G: multivariate analysis of genetic variation in multiple populations." Heredity 112, no. 1 (2014): 21-29.
#' @export
#' @examples 
#' data(dentus)
#' dentus.matrices = dlply(dentus, .(species), function(x) cov(x[-5]))
#' KrzSubspace(dentus.matrices, k = 2)
#' 
#' \dontrun{
#' # The method in Aguirre et al. 2014 can de implemented usign this function as follows:
#' 
#' #Random input data with dimensions traits x traits x populations x MCMCsamples:
#' cov.matrices = aperm(aaply(1:10, 1, function(x) laply(RandomMatrix(6, 40,                    
#'                                                       variance = runif(6,1, 10)), 
#'                            identity)), 
#'                      c(3, 4, 1, 2))
#'    
#' library(magrittr)
#' Hs = alply(cov.matrices, 4, function(x) alply(x, 3)) %>% llply(function(x) KrzSubspace(x, 3)$H)
#' avgH = Reduce("+", Hs)/length(Hs)
#' avgH.vec <- eigen(avgH)$vectors
#' MCMC.H.val = laply(Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))
#' 
#' # confidence intervals for variation in shared subspace directions
#' library(coda)
#' HPDinterval(as.mcmc(MCMC.H.val))    
#' }
KrzSubspace <- function(cov.matrices, k = NULL){
  if (is.null(k))
    k = round(dim(cov.matrices[[1]])[1]/2 - 1)
  AA_T <- llply(cov.matrices, 
                function(x) {
                  A = eigen(x)$vectors[,1:k]
                  A %*% t(A)})
  H = Reduce("+", AA_T)
  eigen_H = eigen(H)
  angles = ldply(AA_T, function(AA_Ti) aaply(eigen_H$vectors[,1:k], 2, 
                                             function(b) acos(sqrt(t(b) %*% AA_Ti %*% b))) * (180/pi)) 
  if(is.null(names(cov.matrices))) angles = cbind(paste0("m", 1:length(cov.matrices)), angles)
  colnames(angles) <- c(".id", paste0("e", 1:length(eigen_H$values[1:k])))
  list(H = H,
       k_eVals_H = eigen_H$values[1:k], 
       k_eVec_H = eigen_H$vectors[,1:k],
       k_angles = angles)
}
