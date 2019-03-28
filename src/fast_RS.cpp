#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat RS(arma::mat x, arma::mat y, int num_vectors) {
  int traits = x.n_cols, k;
  arma::vec comparisons(num_vectors);
  arma::vec null_dist(num_vectors);
  arma::vec base_vector = arma::randn(traits, 1);
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

// [[Rcpp::export]]
arma::vec delta_z_corr(arma::mat x, arma::mat y, int num_vectors, arma::mat random_vectors) {
  int k;
  arma::vec comparisons(num_vectors);
  arma::mat dz_x = x * random_vectors;
  arma::mat dz_y = y * random_vectors;
  for(k = 0; k < num_vectors; ++k){
    comparisons(k) = arma::norm_dot(dz_x.col(k), dz_y.col(k));
  }
  return comparisons;
}
