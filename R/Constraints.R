Constraints <-
function (beta, cov.matrix) return (abs (t (Normalize (eigen (cov.matrix)$vectors[,1])) %*% Normalize (cov.matrix %*% beta)))
