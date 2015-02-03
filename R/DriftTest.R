#'Test drift hypothesis
#'
#'Given a set of covariance matrices and means for terminals, test the hypothesis
#'that obseved divergency is larger/smaller than expected by drift alone.
#'
#'@param means list or array of species means being compared. array must have means in the rows.
#'@param cov.matrix ancestral covariance matrix for all populations
#'@param verbose If TRUE all results are returned, is FALSE only regression coeficient
#'
drift <- function(means, cov.matrix, verbose=TRUE)
{
  if(is.data.frame(means) | (!is.array(means) & !is.list(means)))
    stop("means must be in a list or an array.")
  if(!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if(is.list(means)){
    mean.array = laply(means, identity)
    rownames(mean.array) <- names(means)
  }
  else mean.array <- means
  W.pc <- eigen(cov.matrix)
  projection.Wpc <- as.matrix(mean.array) %*% W.pc$vectors #projecting the means in the principal components of W
  B <- apply(projection.Wpc,2,var) #variance between groups
  log.B <- log(B)
  log.W.pc <- log(W.pc$values) 
  regression <- lm(log.B~log.W.pc)
  regression.coef <- coef(regression)
  CI.95 <- confint(regression)
  CI.regression <- predict(regression, interval="confidence")   
  if(verbose == TRUE){
  objeto <- list(B= log.B, "projecao_medias_W"= projection.Wpc,  
                 "Coeficiente"=regression.coef,"Confidence interval 95%"=CI.95,"PCA_de_W"=W.pc,"Vetor B"=B)
  }
  else{
    objeto <- list("Coeficiente"=regression.coef,"Confidence interval 95%"=CI.95)
  }
  quartz()
  plot(log.B~ log.W.pc, xlab="log(W Eigenvalues)", ylab= "log(B variances)", las=1,
       type= "n")
  text(log.W.pc, log.B,  labels=as.character(1:ncol(cov.matrix)), cex=0.8)  
  abline(regression) #reta de regressao
  lines(log.W.pc, CI.regression[,2], lty=2, col= "red")
  lines(log.W.pc, CI.regression[,3], lty=2, col= "red")
return(objeto)
}

drift(medias.morfologia.reid, cov.matrix= w.morfologia, verbose= FALSE)
drift(medias.morfologia.reid, cov.matrix= w.morfologia, verbose= TRUE)
########################################################
