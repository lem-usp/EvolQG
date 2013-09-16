hh.mod <-
function (mat, hip, nsk = 10000)
{
  with (load.of.functions,
        {
          out <- HansenHouleAverage (mat, nsk)
          HansenHouleWrap2 <- function (hh.func)
            return (apply (hip, 2, hh.func, cov.matrix = mat))
          hip <- apply (hip, 2, Normalize)
          out$mod <- sapply (hansen.houle, HansenHouleWrap2)
          return (out)
        })
}
