#' Print Matrix to file
#'
#' Print a matrix or a list of matrices to file
#' @param  x Matrix or list of matrices
#' @param  out.file Output file
#' @param  ... Aditional parameters
#' @return Prints coma separated matrices, with labels
#' @export
#' @rdname PrintMatrix
#' @author Diogo Melo
#' @examples
#' m.list <- RandomMatrix(10, 4)
#' PrintMatrix(m.list)
PrintMatrix <- function(x, ...) UseMethod('PrintMatrix')

#' @rdname PrintMatrix
#' @method PrintMatrix default
#' @S3method PrintMatrix default
PrintMatrix.default <- function(x, out.file = './matrix.csv', ...){
    write.csv(x, out.file)
}

#' @rdname PrintMatrix
#' @method PrintMatrix list
#' @S3method PrintMatrix list
PrintMatrix.list <- function(x, out.file = './matrix.csv', ...){
    if(is.null(names(x))) names(x) <- 1:length(x)
    sink(out.file, type="output")
    invisible(
              lapply(names(x), function(mat) {
                     y <- x[names(x) == mat]
                     cat(mat)
                     cat('\n')
                     dump(write.table(y, row.names = FALSE, col.names = FALSE, sep=","))
                     cat('\n')
}))
    sink()
}
