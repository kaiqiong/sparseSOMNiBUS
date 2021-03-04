#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "utils.h"



NumericVector proximalOperatorCpp(const double& t, 
                                  const double& lambda1, 
                                  const arma::colvec& u_p,
                                  const int& nk);

List thetaUpdateCpp(const double& stepSize,
                    const NumericVector& theta, 
                    const NumericVector& gBinomLoss, 
                    const int& nk,
                    const int& numCovs, 
                    const double& lambda1, 
                    const arma::mat& basisMat0, 
                    const DataFrame& dat, 
                    const List& designMat1,
                    const NumericVector& theta_m,
                    const int& iter,
                    const bool& accelrt,
                    const bool& truncation);

List oneUpdateCpp(const NumericVector& theta,
               double stepSize, 
               const double& lambda1, 
               const DataFrame& dat,
               const arma::mat& basisMat0,
               const int& nk,
               const int& numCovs,
               const double& shrinkScale,
               const List& designMat1, 
               const NumericVector& theta_m, 
               const int& iter, 
               const bool& accelrt,
               const bool& truncation);


