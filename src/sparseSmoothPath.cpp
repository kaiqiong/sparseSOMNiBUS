#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "utils.h"
#include "updates.h"
#include "proxGradFit.h"

// [[Rcpp::export]]

// run the proximal gradient descent iteration for a sequence of lambda1 given lambda 2

List fitProxGradCppSeq(const NumericVector& ulam,
                    NumericVector& theta, 
                    double& intStepSize,
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
                    const bool& truncation,
                    const double& neg2loglikSat
                    ){
  
  int nlam = ulam.length();
  int myp = (numCovs+1)*nk;
  
  arma::mat thetaMat(nlam, myp);
  arma::mat gNeg2loglikTilda(nlam, myp);

  double lamnow = ulam[0]+100000;
  
  List fit1 = fitProxGradCpp(theta, intStepSize,lamnow, dat, basisMat0, nk,Hp,
                         maxInt, epsilon, shrinkScale,
                         accelrt, numCovs, designMat1,  truncation);
  arma::rowvec temp1 = fit1["thetaEst"];
  thetaMat.row(0) = temp1;
  
  
  arma::rowvec temp2 = fit1["gNeg2loglik"];
  gNeg2loglikTilda.row(0) = temp2;
  // check the fit --- can be removed later -- check on R instead
  
 // NumericVector temp2 = optimcheck(Named("ssfit", fit1), 
  //                                Named("lambda1",lamnow),Named("Hp", Hp), Named("L", L), Named("Hpinv", Hinv),
  //                                Named("nk",nk), );
 // arma::rowvec temp3 = temp2;
 // checkall.row(0) = temp3;
  // end of check
  for (int i=1; i < nlam; ++i){
    NumericVector newInt = fit1["thetaEst"];
    fit1 = fitProxGradCpp(newInt, intStepSize,ulam[i], dat, basisMat0, nk,Hp,
                           maxInt, epsilon, shrinkScale,
                           accelrt, numCovs, designMat1,  truncation);
    arma::rowvec temp1 = fit1["thetaEst"];
    thetaMat.row(i) = temp1;
    
    arma::rowvec temp2 = fit1["gNeg2loglik"];
    gNeg2loglikTilda.row(i) = temp2;
    
    double val1 = fit1["neg2loglik"];
    
    if(val1 < neg2loglikSat) break;
      
  }
  List output=List::create(Named("thetaMat")=thetaMat, 
                           Named("gNeg2loglikTilda") = gNeg2loglikTilda);
  return(output);

}

