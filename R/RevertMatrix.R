#' Revert Matrix
#' 
#' Constructs a covariance matrix based on scores over covariance matrix eigentensors.
#' 
#' @param values vector of values to build matrix, each value corresponding 
#' to a score on the ordered set of eigentensors (up to the maximum number of 
#' eigentensors on the target decomposition)
#' @param etd Eigentensor decomposition of m covariance matrices for 
#' k traits (obtained from \code{\link{EigenTensorDecomposition}})
#' @param scaled should we treat each score as a value given in standard 
#' deviations for each eigentensor? Defaults to TRUE
#' 
#' @return A symmetric covariance matrix with k traits
#' 
#' @references Basser P. J., Pajevic S. 2007. Spectral decomposition of a 4th-order 
#' covariance tensor: Applications to diffusion tensor MRI. Signal Processing. 87:220–236.
#' @references Hine E., Chenoweth S. F., Rundle H. D., Blows M. W. 2009. Characterizing 
#' the evolution of genetic variance using genetic covariance tensors. Philosophical 
#' transactions of the Royal Society of London. Series B, Biological sciences. 364:1567–78.
#' 
#' @examples 
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
      values <- values * etd $ values
    
    rebuilt <- aaply (1:n.dim, 1, function (i)
      etd $ matrices [, , i] * values [i])
    
    rebuilt <- aaply (rebuilt, c(2, 3), sum)
    
    etd $ mean.sqrt %*% expm (rebuilt) %*% etd $ mean.sqrt
  }