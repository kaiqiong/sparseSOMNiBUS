#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "utilsClean.h"




NumericVector proximalOperatorCpp(const double& t, 
                                  const double& lambda1, 
                                  const arma::colvec& u_p,
                                  const int& nk);



double dot_std(NumericVector& x,
               NumericVector& y);



List thetaUpdateCpp(const double& stepSize,
                    const NumericVector& theta, 
                    const NumericVector& gBinomLoss, 
                    const int& nk,
                    const int& numCovs, 
                    const double& lambda1, 
                    const NumericVector& theta_m,
                    const int& iter,
                    const bool& accelrt);


List oneUpdateCpp(const NumericVector& theta,
                  const double& current_lossval,
                  const NumericVector& gBinomLossNum,
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


List binomObjectLossOnlyCpp(const NumericVector& theta,
                            const arma::mat& basisMat0, 
                            const DataFrame& dat, 
                            const int& nk,
                            const int& numCovs,
                            const List& designMat1, 
                            const bool& truncation);