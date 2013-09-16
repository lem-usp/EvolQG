Flexibility <-
function (beta, cov.matrix) return (t (beta) %*% cov.matrix %*% beta / Norm (cov.matrix %*% beta))
