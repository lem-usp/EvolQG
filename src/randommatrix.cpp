#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using Eigen::MatrixXd;

// [[Rcpp::export]]
Eigen::MatrixXd createRandomMatrix(int dimension, float ke) {
    MatrixXd randomMatrix = MatrixXd::Zero(dimension, dimension) ;
    MatrixXd b = MatrixXd::Zero(dimension, dimension);
    b.triangularView<Eigen::Lower>().fill(1.0);
    //MatrixXd fixed = MatrixXd::Ones(dimension, dimension); //what is this line?
    
    randomMatrix.block(1, 0, dimension - 1, 1) = Eigen::MatrixXd::Random(dimension - 1, 1); //contains call to the 
                                                                                            // C random number generator;
                                                                                            //must be replaced by an R
                                                                                            //equivalent for CRAN
    b.block(1, 0, dimension - 1, 1) = randomMatrix.block(1, 0, dimension - 1, 1);

    for (int i = 1; i < dimension; ++i)
      b.block(i, 1, 1, i).fill(sqrt(1 - pow(randomMatrix(i, 0), 2)));

    for (int i = 2; i < dimension; ++i) {
      for (int j = 1; j < i; ++j) {
        float b1, b2, z, y, cosinv, sinTheta;
        Eigen::VectorXd b1a = b.row(j).segment(0, j);
        Eigen::VectorXd b1b = b.row(i).segment(0, j);

        b1 = b1a.dot(b1b);
        b2 = b(j, j)*b(i, j);
        z = b1 + b2;
        y = b1 - b2 ;
        if (b2 < ke){
          randomMatrix(i, j) = b1;
        } else {
          randomMatrix(i, j) = y + (z - y)*R::runif(0, 1);
        }
        
        cosinv = (randomMatrix(i, j) - b1)/b2;
        
        if (b2 != 0) {
          if (cosinv > 1) {
            b.block(i, j + 1, 1, dimension - j - 2).fill(0);
          } else if (cosinv < -1) {
            b(i, j) = -b(i, j);
            b.block(i, j + 1, 1, dimension - j - 2).fill(0);
          } else {
            b(i, j) = b(i, j)*cosinv;
            sinTheta = sqrt(1 - pow(cosinv, 2));
            b.block(i, j + 1, 1, dimension - j - 2) *= sinTheta;
          }
        } 
      }
      
    }
    
    randomMatrix = randomMatrix + randomMatrix.transpose() + MatrixXd::Identity(dimension, dimension);
    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(dimension);
    perm.setIdentity();
    randomMatrix = perm.transpose() * randomMatrix * perm;

    return randomMatrix;
}
