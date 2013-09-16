MeanSquaredCorrelation <-
function (cov.matrix) return (mean (cov2cor (cov.matrix) [lower.tri (diag (nrow (cov.matrix)))]^2))
