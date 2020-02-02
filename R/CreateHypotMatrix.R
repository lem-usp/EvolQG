#' Creates binary correlation matrices
#' 
#' Takes a binary vector or column matrix and generates list of binary correlation matrices representing
#' the partition in the vectors. 
#' @param modularity.hypot Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @return binary matrix or list of binary matrices. If a matrix is passed, all the vectors are combined in the 
#' last binary matrix (total hypothesis of full integration hypotesis).
#' @export
#' @examples
#' rand.hypots <- matrix(sample(c(1, 0), 30, replace=TRUE), 10, 3)
#' CreateHypotMatrix(rand.hypots) 
CreateHypotMatrix <- function(modularity.hypot){
  if(is.null(dim(modularity.hypot))) return(outer(modularity.hypot, modularity.hypot))
  num.hyp <- dim (modularity.hypot) [2]
  num.traits <- dim (modularity.hypot) [1]
  m.hyp.list <- alply(modularity.hypot, 2, function(x) outer(x, x))
  m.hyp.list[[num.hyp+1]] <- matrix(as.integer (as.logical (Reduce ("+", m.hyp.list[1:num.hyp]))),
                                    num.traits, num.traits, byrow=T)
  return(m.hyp.list[1:(num.hyp+1)])
}

CombineHypot <- function(modularity.hypot){
  n.hypots = dim(modularity.hypot)[2]  
  if(is.null(n.hypots)) { # if single hypothesis
    cor.hypot = CreateHypotMatrix(modularity.hypot)
    diag(cor.hypot) <- 1
    return(cor.hypot)
  }
  if(is.null(colnames(modularity.hypot))) colnames(modularity.hypot) <- 1:n.hypots
  counter = BinToDec(rep(1, n.hypots))
  hypot_list = list(null = diag(dim(modularity.hypot)[1]))
  k = 2
  for(i in seq(counter)){
    mask = DecToBin(i)
    mask = as.logical(as.numeric((mask[(32-(n.hypots-1)):32])))
    if(sum(mask) > 1) new_hypot = CreateHypotMatrix(modularity.hypot[,mask])[[sum(mask)+1]]
    else new_hypot = CreateHypotMatrix(modularity.hypot[,mask])
    diag(new_hypot) <- 1
    if(!any(laply(hypot_list, function(x) all(x == new_hypot)))){ 
      hypot_list[[k]] = new_hypot
      names(hypot_list)[[k]] <- paste(colnames(modularity.hypot)[mask], collapse = "_")
      k = k + 1
    }
  }
  hypot_list
}

#' Create binary hypothesis
#' 
#' Takes a vetor describing a trait partition and returns a binary matrix of the partitions where each line represents a trait and each column a module. In the output matrix, if modularity.hypot[i,j] == 1, trait i is in module j.
#' @param x vector of trait partition. Each partition receive the same symbol.
#' @return Matrix of hypothesis. Each line represents a trait and each column a module.
#' if modularity.hypot[i,j] == 1, trait i is in module j.
#' @export
#' @examples
#' x = sample(c(1, 2, 3), 10, replace = TRUE)
#' Partition2HypotMatrix(x) 
Partition2HypotMatrix <- function(x){
  sapply(unique(x), function(i) as.numeric(x == i))  
}

# http://stackoverflow.com/questions/12892348/convert-binary-string-to-binary-or-decimal-value
BinToDec <- function(x) sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
# http://stackoverflow.com/questions/6614283/converting-decimal-to-binary-in-r
DecToBin <- function(x) sapply(strsplit(paste(rev(intToBits(x))),""),`[[`,2)