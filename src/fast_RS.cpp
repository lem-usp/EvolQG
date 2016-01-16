#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;

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
NumericVector normalize(NumericVector x) {
  int n = x.size();
  x = x / sqrt(sum(pow((x), 2)));
  return x;
}

// [[Rcpp::export]]
NumericVector RS(NumericMatrix x, NumericMatrix y, int num_vectors) {
  int traits = x.nrow();
  int i, j, k;

  NumericVector out(3);
  NumericVector base_vector(traits);
  NumericVector random_vector(traits);
  NumericVector dz_x(traits);
  NumericVector dz_y(traits);
  NumericVector null_dist(num_vectors);
  NumericVector comparisons(num_vectors);
  base_vector = normalize(rnorm(traits));
  
  for(k = 0; k < num_vectors; ++k){
    random_vector = normalize(rnorm(traits));
    null_dist[k] = sum(random_vector * base_vector);
    for(i = 0; i < traits; ++i){
      dz_x[i] = 0;
      dz_y[i] = 0;
      for(j = 0; j < traits; ++j){
        dz_x[i] += x(i, j) * random_vector[j];
        dz_y[i] += y(i, j) * random_vector[j];
      }
      comparisons[k] = sum(normalize(dz_x) * normalize(dz_y));
    }
  }
  
  out[0] = mean(comparisons);
  out[1] = 0;
  for(k = 0; k < num_vectors; k++){
    if(null_dist[k] > out[0]) out[1]++;
  }
  out[1] /= num_vectors;
  out[2] = sd(comparisons);
  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(evolqg)
library(microbenchmark)
x = RandomMatrix(40)
y = RandomMatrix(40)
microbenchmark(
  rscpp = RS(x, y, 1000),
  rsR = RandomSkewers(x, y)
)
*/
