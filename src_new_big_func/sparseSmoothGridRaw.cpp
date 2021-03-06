#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "utilsClean.h"
#include "updatesClean.h"
#include "proxGradFitClean.h"




//'@title fitProxGradGridRaw.cpp
//'@description given a grid of lambda and a sequence of alpha fit the sequence for each alpha
//'@param lamGrid A grid values of lambda1, each column corresponds to a different ulam2 
//'@param ulam2 a sequence of alpha/lambda2
//'@param Linv a list for the matrix Linv; its length equal to length(ulam2); inverse of the cholesky decomp L; H = t(L)%*%L
//'@param designMat1 design matrix for the original version


// [[Rcpp::export]]


List sparseSmoothGridRawCpp(const arma::mat& lamGrid,
                            const NumericVector& ulam2,
                            const List& Linv,
                            NumericVector& theta, 
                            double& intStepSize,
                            const DataFrame& dat, 
                            const arma::mat& basisMat0, 
                            const int& nk, 
                            const double& maxInt,
                            const double& epsilon,
                            const double& shrinkScale,
                            const bool& accelrt, 
                            const int& numCovs, 
                            const List& designMat1, 
                            const bool& truncation,
                            const double& neg2loglikSat){
  
  int nlam2 = ulam2.length();
  int myp = (numCovs+1)*nk;
  int nlam = lamGrid.n_rows;
  
  // Initialization:
  
  arma::cube thetaMat(myp, nlam, nlam2);
  arma::cube thetaMatOriginal(myp, nlam, nlam2);
  
  
  for (int ii=0; ii < nlam2; ++ii){ // loop over different ulam2, different columns of lamGrid
    
    arma::vec ulam = lamGrid.col(ii); // the sequence of lambda1 for warm start
    
    
    // step 1 calculate design tilda
    
    arma::mat tt = Linv[ii];
    
    arma::mat basisMat0_tilda = basisMat0 * tt;
    
    List designMat1_tilda(numCovs);
    for (int jj =0; jj <numCovs; ++jj){
      arma::mat temp = designMat1[jj];
      designMat1_tilda[jj] = temp * tt ;
    }
    
    
    // step 2 fit a sequence of ulam
    
    List fit1 = fitProxGradCppClean(theta, intStepSize, lamGrid(0,ii), dat, basisMat0_tilda, nk,
                                    maxInt, epsilon, shrinkScale,
                                    accelrt, numCovs, designMat1_tilda,  truncation);
    arma::vec temp1 = fit1["thetaEst"];  // estimate of theta_tilda
    
    arma::vec thetaOri(myp); // the estimate of the original theta
    
    List temp1sep = fit1["thetaEstSep"];
    
    arma::vec nowsep = temp1sep[0];
    thetaOri(span(0, nk-1)) = tt * nowsep;
    
    for(int ll=0; ll<numCovs; ++ll){
      arma::vec nowsep = temp1sep[ll+1];
      
      thetaOri(span((ll+1)*nk,(ll+2)*nk-1 )) =tt * nowsep;
    }
    
    
    
    thetaMat.subcube(0,0,ii, myp-1,0,ii) = temp1;
    thetaMatOriginal.subcube(0,0,ii, myp-1,0,ii) = thetaOri;
    
    
    for (int i=1; i < nlam; ++i){
      NumericVector newInt = fit1["thetaEst"];
      fit1 = fitProxGradCppClean(newInt, intStepSize, lamGrid(i,ii), dat, basisMat0_tilda, nk,
                                 maxInt, epsilon, shrinkScale,
                                 accelrt, numCovs, designMat1_tilda,  truncation);
      
      arma::vec temp1 = newInt;  // estimate of theta_tilda
    //  arma::vec thetaOri(myp); // the estimate of the original theta
      
      List temp1sep = fit1["thetaEstSep"];
      
      arma::vec nowsep = temp1sep[0];
      thetaOri(span(0, nk-1)) = tt * nowsep;
      
      for(int ll=0; ll<numCovs; ++ll){
        arma::vec nowsep = temp1sep[ll+1];
        
        thetaOri(span((ll+1)*nk,(ll+2)*nk-1 )) =tt * nowsep;
      }
      
      
      
      
      thetaMat.subcube(0,i,ii, myp-1,i,ii) = temp1;
      thetaMatOriginal.subcube(0,i,ii, myp-1,i,ii) = thetaOri;
      
      double val1 = fit1["neg2loglik"];
      
      if(val1 < neg2loglikSat) break;
      
    }
    
    
    
  }
  
  
  
  List output=List::create(Named("thetaMat")=thetaMat, 
                           Named("thetaMatOriginal") = thetaMatOriginal);
  return(output);
  
  
}

