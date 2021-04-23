#' EvolQG
#'
#' The package for evolutionary quantitative genetics.
#'
#' @name evolqg
#' @docType package
#' @import plyr
#' @importFrom graphics abline arrows axis layout mtext par plot text
#' @importFrom grDevices rgb
#' @importFrom methods is
#' @importFrom stats anova confint cor cor.test cov cov2cor lm mahalanobis princomp quantile reorder residuals rnorm runif sd var
#' @importFrom utils write.csv write.table
NULL

#' Linear distances for five mouse lines
#'
#' Skull distances measured from landmarks in 5 mice lines: 4 body weight selection lines and 1 control line. 
#' Originally published in Penna, A., Melo, D. et. al (2017) 10.1111/evo.13304
#'
#' @name ratones
#'
#' @docType data
#'
#' @usage data(ratones)
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references Penna, A., Melo, D., Bernardi, S., Oyarzabal, M.I. and Marroig, G. (2017), The evolution of phenotypic integration: How directional selection reshapes covariation in mice. Evolution, 71: 2370-2380. https://doi.org/10.1111/evo.13304
#' (\href{https://pubmed.ncbi.nlm.nih.gov/28685813/}{PubMed})
#'
#' @source \href{https://datadryad.org/stash/dataset/doi:10.5061/dryad.5gr8r}{Dryad Archive}
#'
#' @examples
#' data(ratones)
#'    
#' # Estimating a W matrix, controlling for line and sex
#' model_formula = paste0("cbind(", 
#'                        paste(names(ratones)[13:47], collapse = ", "),
#'                        ") ~ SEX + LIN")
#' ratones_W_model = lm(model_formula, data = ratones)
#' W_matrix = CalculateMatrix(ratones_W_model)
#' 
#' # Estimating the divergence between the two direction of selection
#' delta_Z = colMeans(ratones[ratones$selection == "upwards", 13:47]) -
#'           colMeans(ratones[ratones$selection == "downwards", 13:47])
#'           
#'  # Reconstructing selection gradients with and without noise control         
#' Beta = solve(W_matrix, delta_Z)
#' Beta_non_noise = solve(ExtendMatrix(W_matrix, ret.dim = 10)$ExtMat, delta_Z)
#' 
#' # Comparing the selection gradients to the observed divergence
#' Beta %*% delta_Z /(Norm(Beta) * Norm(delta_Z))
#' Beta_non_noise %*% delta_Z /(Norm(Beta_non_noise) * Norm(delta_Z))      
#'           
"ratones"

#' Example multivariate data set
#'
#' Simulated example of 4 continuous bone lengths from 5 species.
#'
#' \itemize{
#' \item humerus 
#' \item ulna 
#' \item femur 
#' \item tibia 
#' \item species 
#' }
#'
#' @docType data
#' @keywords datasets
#' @name dentus
#' @usage data(dentus)
#' @format A data frame with 300 rows and 5 variables
NULL

#' Tree for dentus example species
#'
#' Hypothetical tree for the species in the dentus data set.
#'
#' @docType data
#' @keywords datasets
#' @name dentus.tree
#' @usage data(dentus.tree)
#' @format ape tree object
NULL