#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "utils.h"
#include "updates.h"








// [[Rcpp::export]]

List fitProxGradCpp(NumericVector& theta, 
                    double& intStepSize,
                    const double& lambda1,
                    const DataFrame& dat, 
                    const arma::mat& basisMat0, 
                    const int& nk, 
                    const double& maxInt,
                    const double& epsilon,
                    const double& shrinkScale,
                    const bool& accelrt, 
                    const int& numCovs, 
                    const List& designMat1, 
                    const bool& truncation){
  // run the proximal gradient descent iteration for a given pairs of lambda1 and lambda2.
  
  // arma::mat Hp=sparOmega + lambda2*smoOmega1; 
  
 // arma::mat Hp=(1-lambda2)*sparOmega + lambda2*smoOmega1; 
  int iter=0;
  //-- for --test//
  
  int nsamp = dat.nrows();

 // arma::mat thetaMat(maxInt+10, theta.length());
  
  NumericVector stepSizeVec(maxInt+1);
  //-- for --test//
  List out=oneUpdateCpp(theta, intStepSize, lambda1, dat, basisMat0, nk,   numCovs, shrinkScale, 
                        designMat1,theta, iter, accelrt, truncation);
  
  NumericVector theta_new = out["theta_l_proximal"];
  stepSizeVec(iter) = out["stepSize"];

  
  arma::vec diff=theta_new-theta;
  double innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double tol = pow(innerdot, 0.5);
  
  // Calculate the loss value evaluated at the theta_new
  List binom_out_new = out["binom_out_new"];
  double neg2loglik = binom_out_new["neg2loglik"];
  List theta_l_proximal_sep=out["theta_l_proximal_sep"];
  double penTerms = twoPenaltiesCpp(theta_l_proximal_sep, lambda1, numCovs, nk);
  double lossSum = neg2loglik + penTerms;
  double lossSumOld = lossSum;
  // 
  
  double tol2 = 100;
  
 // arma::rowvec temp = theta_new;
 // thetaMat.row(iter) = temp;

  
  while(((tol > epsilon) & (iter < maxInt-1)) & (tol2 >epsilon)){
  
   lossSumOld = lossSum;
   theta = theta_new;
 
   out=oneUpdateCpp(theta, intStepSize, lambda1, dat, basisMat0, nk,  numCovs, shrinkScale, 
                     designMat1,theta, iter, accelrt, truncation);
   theta_new=out["theta_l_proximal"];
 
   
   // Rcout << "step size : " << intStepSize << "\n";
   
    diff=theta_new-theta;
    innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    tol = pow(innerdot, 0.5);
    
  // Add another tolerance checking for the values of objective functions  
     theta_l_proximal_sep=out["theta_l_proximal_sep"];
     binom_out_new = out["binom_out_new"];
     neg2loglik = binom_out_new["neg2loglik"];
    
     penTerms = twoPenaltiesCpp(theta_l_proximal_sep, lambda1, numCovs, nk);
     lossSum = neg2loglik + penTerms;
    
    tol2 = lossSumOld-lossSum;
    tol2 = tol2/nsamp;
  // 
    iter = iter+1;
    stepSizeVec(iter) = out["stepSize"];
   // arma::rowvec temp = theta_new;
   // thetaMat.row(iter) = temp;
   

  }
 
// Rcout << "Iter: " << iter << "\n";
 
 // do{
 //  List out=oneUpdateCpp(theta, stepSize, lambda1, dat, basisMat0, nk, Hp,  numCovs, shrinkScale, 
 //                         designMat1,theta, iter, accelrt, truncation);
  
  //  NumericVector theta_new = out["theta_l_proximal"];
  //  NumericVector diff=theta_new-theta;
  //  double innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  //  double tol = pow(innerdot, 0.5);
  //  iter = iter+1;
  //  theta = theta_new;
 // } while(tol > epsilon & iter < maxInt);
  
  //NumericVector gNeg2loglik= binom_out_new["gNeg2loglik"];
  
 
                                                                        
  List output=List::create(Named("thetaEst")=theta_new, 
                        Named("lossSum")=lossSum,
                        Named("neg2loglik")=neg2loglik,
                        Named("thetaEstSep") = theta_l_proximal_sep,
                      //  Named("pi_ij")=pi_ij,
                        Named("Iter") = iter,
                    //    Named("gNeg2loglik")=gNeg2loglik,
                        Named("stepSizeVec") = stepSizeVec);
  return(output);
        
}


// [[Rcpp::export]]
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
                    const bool& truncation){
  
   arma::mat Hp=sparOmega + lambda2*smoOmega1; 
  
 // arma::mat Hp=(1-lambda2)*sparOmega + lambda2*smoOmega1; 
  int iter=0;
  //-- for --test//
  
  // arma::mat thetaMat(maxInt+10, theta.length());
  
  //-- for --test//
  List out=oneUpdateCpp(theta, intStepSize, lambda1, dat, basisMat0, nk,  numCovs, shrinkScale, 
                        designMat1,theta, iter, accelrt, truncation);
  
  NumericVector theta_new = out["theta_l_proximal"];
  
  
  
  arma::vec diff=theta_new-theta;
  double innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double tol = pow(innerdot, 0.5);
  
  // arma::rowvec temp = theta_new;
  // thetaMat.row(iter) = temp;
  
  while((tol > epsilon) & (iter < maxInt)){
    
    theta = theta_new;
    
    
    
    out=oneUpdateCpp(theta, intStepSize, lambda1, dat, basisMat0, nk,   numCovs, shrinkScale, 
                     designMat1,theta, iter, accelrt, truncation);
    theta_new=out["theta_l_proximal"];
    
    
    //  Rcout << "step size : " << intStepSize << "\n";
    
    diff=theta_new-theta;
    innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    tol = pow(innerdot, 0.5);
    iter = iter+1;
    
    // arma::rowvec temp = theta_new;
    // thetaMat.row(iter) = temp;
    
    
  }
  
  // Rcout << "Iter: " << iter << "\n";
  
  // do{
  //  List out=oneUpdateCpp(theta, stepSize, lambda1, dat, basisMat0, nk, Hp,  numCovs, shrinkScale, 
  //                         designMat1,theta, iter, accelrt, truncation);
  
  //  NumericVector theta_new = out["theta_l_proximal"];
  //  NumericVector diff=theta_new-theta;
  //  double innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  //  double tol = pow(innerdot, 0.5);
  //  iter = iter+1;
  //  theta = theta_new;
  // } while(tol > epsilon & iter < maxInt);
  
  
  
  List theta_l_proximal_sep=out["theta_l_proximal_sep"];
  double penTerms = twoPenaltiesCpp(theta_l_proximal_sep,  lambda1, numCovs, nk);
  List binom_out_new = out["binom_out_new"];
  NumericVector pi_ij = binom_out_new["pi_ij"];
  double neg2loglik = binom_out_new["neg2loglik"];
  double lossSum = neg2loglik + penTerms;
  
  List thetaEstSep = out["theta_l_proximal_sep"];                                                                           
  List output=List::create(Named("thetaEst")=theta_new, 
                           Named("lossSum")=lossSum,
                           Named("thetaEstSep") = thetaEstSep,
                           Named("pi_ij")=pi_ij);
  return(output);
  
}

