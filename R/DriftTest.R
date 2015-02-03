#'Test drift hypothesis
#'
#'Given a set of covariance matrices and means for terminals, test the hypothesis
#'that obseved divergency is larger/smaller than expected by drift alone.
#'
#'@param means list or array of species means being compared. array must have means in the rows.
#'@param cov.matrix ancestral covariance matrix for all populations
#'@param show.plot boolean. If TRUE, plot of eigenvalues of ancetral matrix by between group variance is showed.
#'@return list of results containing:
#'@return regresion: the linear regression between the log of the eigenvalues of the ancestral matrix and the log of the between group variance (projected on the eigenvectors of the ancenstral matrix)
#'@return coefficient_CI_95: confidence intervals for the regression coefficients
#'@return log.between_group_variance: log of the between group variance (projected on the ancestral matrix eigenvectors)
#'@return log.W_eVals: log of the ancestral matrix eigenvalues
#'@return plot: plot of the regression using ggplot2
#'@note If the regression coefficient is significantly different to one, the null hypothesis of drift is rejected.
#'@references Marroig, G., and Cheverud, J. M. (2004). Did natural selection or genetic drift 
#'produce the cranial diversification of neotropical monkeys? The American Naturalist, 163(3), 417-428. doi:10.1086/381693
#'@references Prôa, M., O'Higgins, P. and Monteiro, L. R. (2013), Type I error rates for testing genetic drift with phenotypic covariance matrices: A simulation study. Evolution, 67: 185-195. doi: 10.1111/j.1558-5646.2012.01746.x
#'@author Ana Paula Assis, Diogo Melo
#'@export
#'@import plyr
#'@importFrom ggplot2 ggplot geom_text geom_smooth labs theme_bw
#'@examples
#'
#'means = array(rnorm(40*10), c(10, 40)) 
#'cov.matrix = RandomMatrix(40, 1, 1, 10)
#'DriftTest(means, cov.matrix)
DriftTest <- function(means, cov.matrix, show.plot=TRUE)
{
  if(is.data.frame(means) | (!is.array(means) & !is.list(means)))
    stop("means must be in a list or an array.")
  if(!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if(is.list(means)){
    mean.array <- laply(means, identity)
  }  else{ 
    mean.array <- means
  }
  W.pc <- eigen(cov.matrix)
  projection.Wpc <- as.matrix(mean.array) %*% W.pc$vectors #projecting the means in 
                                                           #the principal components of W
  log.B_variance <- log(apply(projection.Wpc, 2, var)) #variance between groups
  log.W_eVals <- log(W.pc$values) 
  regression <- lm(log.B_variance~log.W_eVals)
  reg.plot <- ggplot(data.frame(log.B_variance, 
                                log.W_eVals, 
                                names = 1:(dim(mean.array)[2])), 
                     aes(log.W_eVals, log.B_variance)) + 
              geom_text(aes(label = names)) + 
              geom_smooth(method = "lm", color = 'black') + 
              labs(x = "log(W Eigenvalues)", y = "log(B variances)") + theme_bw()
  if(show.plot) print(reg.plot)
  objeto <- list("regression" = regression,
                 "coefficient_CI_95" = confint(regression),
                 "log.between_group_variance" = log.B_variance, 
                 "log.W_eVals" = log.W_eVals,
                 "plot" = reg.plot)
  return(objeto)
}