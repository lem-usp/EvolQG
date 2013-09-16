CalcR2 <-
function (c.matrix){
    cor.matrix = cov2cor(c.matrix)
    return (mean (cor.matrix [lower.tri (cor.matrix)]^2))
}
