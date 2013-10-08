PrintMatrix <- function(x, ...) UseMethod('PrintMatrix')

PrintMatrix.default <- function(x, out.file = './matrix.csv', ...){
    write.csv(x, out.file)
}

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
