Autonomy <-
function (beta, cov.matrix) return (((t (beta) %*% solve (cov.matrix) %*% beta)^(-1)) / (t (beta) %*% cov.matrix %*% beta))
