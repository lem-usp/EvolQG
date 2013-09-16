Pc1Percent <-
function (cov.matrix) return (eigen (cov.matrix)$values [1] / sum (eigen (cov.matrix)$values))
