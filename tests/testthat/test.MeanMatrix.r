test_that("MeanMatrix works on dentus data",
          {
            data(dentus)
            dentus.vcv <- daply (dentus, .(species), function(x) cov(x[,-5]))
            dentus.vcv <- aperm(dentus.vcv, c(2, 3, 1))
            dentus_mean = MeanMatrix(dentus.vcv)
            expect_is(MeanMatrix(dentus.vcv), "matrix")
            expect_is(MeanMatrix(dentus.vcv[,,1:3]), "matrix")
          }
)

test_that("MeanMatrix works on carray data",{
  skip_on_cran()
  carray_data = array(NA, dim = c(5, 5, 3))
  carray_data[,,1] = c(0.8467  ,0.5195  ,0.8089  ,0.2300, -0.0228,
                       0.5195  ,0.4611  ,0.5595  ,0.4020, -0.2681,
                       0.8089  ,0.5595  ,0.9253,  0.5862, -0.2901,
                       0.2300  ,0.4020  ,0.5862  ,1.0472, -0.7404,
                       -0.0228, -0.2681, -0.2901, -0.7404,  2.2454)
  
  carray_data[,,2] = c(0.7395, 0.6315, 0.8940, 0.4118, 0.1810,
                       0.6315, 0.8664, 0.7663, 0.4747, 0.2356,
                       0.8940, 0.7663, 1.1668, 0.7014, 0.2386,
                       0.4118, 0.4747, 0.7014, 0.7541, 0.1393,
                       0.1810, 0.2356, 0.2386, 0.1393, 0.9006)
  
  carray_data[,,3] = c(0.9727 ,0.7276  ,0.9256, 0.2217, -0.0499,
                       0.7276 ,0.8336  ,0.7469, 0.5929,  0.3071,
                       0.9256 ,0.7469  ,1.0325, 0.5313, -0.0542,
                       0.2217 ,0.5929  ,0.5313, 1.0943,  0.4322,
                       -0.0499, 0.3071, -0.0542, 0.4322,  0.9342)
  expect_is(MeanMatrix(carray_data), "matrix")
})
