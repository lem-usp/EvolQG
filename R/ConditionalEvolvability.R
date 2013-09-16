ConditionalEvolvability <-
function (beta, cov.matrix) return ((t (beta) %*% solve (cov.matrix) %*% beta)^(-1))
