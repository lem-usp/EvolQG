#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using Eigen::MatrixXd;

// [[Rcpp::export]]
Eigen::MatrixXd createRandomMatrix(int dimension, float ke) {
    MatrixXd randomMatrix = MatrixXd::Zero(dimension, dimension) ;
    MatrixXd b = MatrixXd::Zero(dimension, dimension);

    b.triangularView<Eigen::Lower>().fill(1.0);
    randomMatrix.block(1, 0, dimension - 1, 1) = Eigen::MatrixXd::Random(dimension - 1, 1);
    b.block(1, 0, dimension - 1, 1) = randomMatrix.block(1, 0, dimension - 1, 1);
    b.triangularView<Eigen::Lower>() = randomMatrix.triangularView<Eigen::Lower>();

    for (int i = 1; i < dimension; ++i)
      b.block(i, 1, 1, i).fill(sqrt(1 - pow(randomMatrix(i, 0), 2)));

    for (int i = 2; i < dimension; ++i) {
      for (int j = 1; j <= i - 1; ++j) {
        float b1, b2, z, y, cosinv, sinTheta;
        Eigen::VectorXd b1a = Eigen::VectorXd::Map(b.block(j, 0, 1, j-1).data(), j-1);
        Eigen::VectorXd b1b = Eigen::VectorXd::Map(b.block(i, 0, 1, j-1).data(), j-1);
        b1 = b1a.dot(b1b);
        b2 = b(j, j)*b(i, j);
        z = b1 + b2;
        y = b1 - b2 ;
        if (b2 < ke){
          randomMatrix(i, j) = b1;
          cosinv = 0;
        } else {
          randomMatrix(i, j) = y + (z - y)*R::runif(0, 1);
          cosinv = (randomMatrix(i, j) - b1)/b2;
        }
        
        if (cosinv > 1) {
          b.block(i, j + 1, 1, dimension - j - 1).fill(0);
        } else if (cosinv < -1) {
          b(i, j) = -b(i, j);
          b.block(i, j + 1, 1, dimension - j - 1).fill(0);
        } else {
          b(i, j) = b(i, j)*cosinv;
          sinTheta = sqrt(1 - pow(cosinv, 2));
          b.block(i, j + 1, 1, dimension - j - 1) *= sinTheta;
        }
        
      }
      
    }
    
    randomMatrix = randomMatrix + randomMatrix.transpose() + MatrixXd::Identity(dimension, dimension);
    // TODO: Rows permutation
    // perm = sample(1:num.traits)
    // random.corr = (random.corr[perm,])[,perm]
    return randomMatrix;
}
