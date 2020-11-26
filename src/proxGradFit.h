#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "utils.h"
#include "updates.h"



// run the proximal gradient descent iteration for a given pairs of lambda1 and lambda2.
List fitProxGradCpp(NumericVector& theta, 
                    double& intStepSize,
                    const double& lambda1,
                    const DataFrame& dat, 
                    const arma::mat& basisMat0, 
                    const int& nk, 
                    const arma::mat& Hp, 
                    const double& maxInt,
                    const double& epsilon,
                    const double& shrinkScale,
                    const bool& accelrt, 
                    const int& numCovs, 
                    const List& designMat1, 
                    const bool& truncation);



List fitProxGradCppold(NumericVector& theta, 
                    double& intStepSize,
                    const double& lambda1,
                    const DataFrame& dat, 
                    const arma::mat& basisMat0, 
                    const int& nk, 
                    const arma::mat& sparOmega, 
                    const double& lambda2, 
                    const arma::mat& smoOmega1, 
                    const double& maxInt,
                    const double& epsilon,
                    const double& shrinkScale,
                    const bool& accelrt, 
                    const int& numCovs, 
                    const List& designMat1, 
                    const arma::mat& basisMat1,
                    const bool& truncation);

