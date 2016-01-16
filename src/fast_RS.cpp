#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::mat RS(arma::mat x, arma::mat y, int num_vectors) {
  int traits = x.n_cols, k;
  arma::vec comparisons(num_vectors);
  arma::vec null_dist(num_vectors);
  arma::vec base_vector = rnorm(traits);
  arma::vec out(3);
  arma::mat random_vectors = arma::randn(traits, num_vectors);
  arma::mat dz_x = x * random_vectors;
  arma::mat dz_y = y * random_vectors;
  for(k = 0; k < num_vectors; ++k){
    comparisons(k) = arma::norm_dot(dz_x.col(k), dz_y.col(k));
    null_dist(k) = arma::norm_dot(base_vector, random_vectors.col(k));
  }
  out(0) = arma::mean(comparisons);
  out(1) = 0;
  for(k = 0; k < num_vectors; k++){
    if(null_dist(k) > out(0)) out(1)++;
  }
  out(1) /= num_vectors;
  out(2) = arma::stddev(comparisons);
  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(evolqg)
library(microbenchmark)
x = RandomMatrix(100)
y = RandomMatrix(100)
RS(x, y, 1000)
RandomSkewers(x, y)
microbenchmark(
  rsarma = RS(x, y, 1000),
  rsR = RandomSkewers(x, y)
)
*/
