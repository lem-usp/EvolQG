#' Revert Matrix
#' 
#' Constructs a covariance matrix based on scores over covariance matrix eigentensors.
#' 
#' @param values vector of values to build matrix, each value corresponding 
#' to a score on the ordered set of eigentensors (up to the maximum number of 
#' eigentensors on the target decomposition); if there are less values than eigentensors 
#' provided in etd (see below), the function will assume zero as the value 
#' for the score in remaining eigentensors
#' @param etd Eigentensor decomposition of m covariance matrices for 
#' k traits (obtained from \code{\link{EigenTensorDecomposition}})
#' @param scaled should we treat each score as a value given in standard 
#' deviations for each eigentensor? Defaults to TRUE
#' 
#' @return A symmetric covariance matrix with k traits
#' 
#' @references Basser P. J., Pajevic S. 2007. Spectral decomposition of a 4th-order 
#' covariance tensor: Applications to diffusion tensor MRI. Signal Processing. 87:220-236.
#' @references Hine E., Chenoweth S. F., Rundle H. D., Blows M. W. 2009. Characterizing 
#' the evolution of genetic variance using genetic covariance tensors. Philosophical 
#' transactions of the Royal Society of London. Series B, Biological sciences. 364:1567-78.
#' 
#' @examples 
#' 
#' ## we can use RevertMatrix to represent eigentensors using SRD to compare two matrices
#' ## which differ with respect to their projections on a single directions
#' 
#' data(dentus)
#' 
#' dentus.vcv <- daply (dentus, .(species), function(x) cov(x[,-5]))
#'
#' dentus.vcv <- aperm(dentus.vcv, c(2, 3, 1))
#'
#' dentus.etd <- EigenTensorDecomposition(dentus.vcv, TRUE)
#' 
#' ## calling RevertMatrix with a single value will use this value as the score
#' ## on the first eigentensor and use zero as the value of remaining scores
#' 
#' low.et1 <- RevertMatrix(-1.96, dentus.etd, TRUE)
#' upp.et1 <- RevertMatrix(1.96, dentus.etd, TRUE)
#'
#' srd.et1 <- SRD(low.et1, upp.et1)
#'  
#' plot(srd.et1)
#'
#' ## we can also look at the second eigentensor, by providing each call 
#' ## of RevertMatrix with a vector of two values, the first being zero
#' 
#' low.et2 <- RevertMatrix(c(0, -1.96), dentus.etd, TRUE)
#' upp.et2 <- RevertMatrix(c(0, 1.96), dentus.etd, TRUE)
#'
#' srd.et2 <- SRD(low.et2, upp.et2)
#'  
#' plot(srd.et2)
#'
#' @export
#' @importFrom expm expm
#' 
#' @rdname RevertMatrix
#' 
RevertMatrix <-
  function (values, etd, scaled = TRUE)
  { 
    n.scores <- length(values)
    n.dim <- length(etd $ values)
    
    if(n.scores < n.dim)
      values <- c(values, rep(0, n.dim - n.scores))
    
    if(n.scores > n.dim)
    {
      warning('More scores provided than actual eigentensors in etd. \n',
              'Using only the first ', n.dim, ' scores.')
      values <- values [1:n.dim]
    }
    
    if(scaled)
      values <- values * sqrt(etd $ values)
    
    rebuilt <- aaply (1:n.dim, 1, function (i)
      etd $ matrices [, , i] * values [i])
    
    rebuilt <- aaply (rebuilt, c(2, 3), sum)
    
    out.mat <- etd $ mean.sqrt %*% expm (rebuilt) %*% etd $ mean.sqrt
    
    if(!is.null(dimnames(etd $ matrices)))
      dimnames(out.mat) <- dimnames(etd $ matrices) [1:2]
    
    out.mat
  }
