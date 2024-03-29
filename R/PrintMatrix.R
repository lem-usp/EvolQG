#' Print Matrix to file
#'
#' Print a matrix or a list of matrices to file
#' @param  x Matrix or list of matrices
#' @param  output.file Output file
#' @param  ... Additional parameters
#' @return Prints coma separated matrices, with labels
#' @export
#' @rdname PrintMatrix
#' @author Diogo Melo
#' @examples
#' m.list <- RandomMatrix(10, 4)
#' tmp = file.path(tempdir(), "matrix.csv")
#' PrintMatrix(m.list, output.file = tmp )
PrintMatrix <- function(x, ...) UseMethod('PrintMatrix')

#' @rdname PrintMatrix
#' @method PrintMatrix default
#' @export
PrintMatrix.default <- function(x, output.file, ...){
    write.csv(x, output.file)
}

#' @rdname PrintMatrix
#' @method PrintMatrix list
#' @export
PrintMatrix.list <- function(x, output.file, ...){
    if(is.null(names(x))) names(x) <- 1:length(x)
    sink(output.file, type="output")
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
