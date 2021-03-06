#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "utilsClean.h"
#include "updatesClean.h"


// [[Rcpp::export]]




//'@title fitProxGradCpp
//'@description use proximial gradient descent with backtracking line search to minimize
//'a penalized negative binomial likelihood with a group LASSO penalty. This function is
//'for the tilda version of our objective function. So matrix Hp is not needed
//'@name fitProxGradCpp
//'@param theta the initial value for theta tilda
//'@param intStepSize the initial step size used
//'@param lambda1 the penalty parameter lambda in our paper
//'@param dat data frame with two columns named as "Meth_Counts" and "Total_Counts"
//'@param basisMat0 design matrix for the intercept. Its row equals to the row of dat
//'@param nk number of knots used for the covariate --- currently it is the same for all
List fitProxGradCppClean(NumericVector& theta, 
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
//Initialization
  
  int iter=0;
  int nsamp = dat.nrows();
 // NumericVector stepSizeVec(maxInt+1);
  
  //----- getSeparateThetaCpp
  int size=numCovs+1;
  List thetaSep(size);
  for (int i=0; i < size; ++i){
    int lower = (i*nk);
    int upper = (i+1)*nk-1;
    thetaSep[i] = theta[myseq(lower, upper)];
  }
  //---done
  
  //---estimatePijCpp
  int myrows=basisMat0.n_rows;
  
  arma::mat lpSep(myrows,size);
  
  arma::vec a = thetaSep[0];
  lpSep.col(0) =  basisMat0 * a; 
  
  for (int i=1; i < size; ++i){
    
    arma::vec b = thetaSep[i];
    
    arma::mat designNow = designMat1[i-1];
    
    lpSep.col(i) = designNow * b;
    
  }

  arma::vec lp_ij  = arma::sum(lpSep,1);

  arma::vec pi_ij = 1/(1+trunc_exp(-lp_ij));
  
  NumericVector  pi_ij_now =NumericVector(pi_ij.begin(), pi_ij.end());
  
  double eps=getPrec(1) * 10;
  if(truncation==TRUE){
    pi_ij_now[ pi_ij_now > (1 - eps)] = 1 - eps;
    pi_ij_now[ pi_ij_now < eps] = eps;
  }
  
  //---done
 
  //---binomObjectCpp and gbinomObjectCpp
  arma::vec meth=dat["Meth_Counts"];
  arma::vec total=dat["Total_Counts"];
  arma::vec unmeth= total-meth;
  
  arma::vec pi_ij_arma=pi_ij_now;
  arma::vec pi_ij_arma_c=1-pi_ij_now;
  arma::vec logpi=arma::log(pi_ij_arma);
  arma::vec log1mpi=arma::log(pi_ij_arma_c);
  
  arma::vec loglik_ij = meth % logpi + unmeth % log1mpi;
  
  double current_lossval =  (-2)*arma::sum(loglik_ij);
      //calculate gradients
  
  arma::vec gPi_ij = meth - total % pi_ij_arma;
  
  int mycols=nk*(numCovs+1);
  
  arma::vec grad(mycols);
  
  arma::mat temp=basisMat0;
  
  temp.each_col() %= gPi_ij;
  arma::mat g0_now= arma::sum(temp,0);// calculat the column sum
  
  grad(span(0, nk-1)) = trans(g0_now);
  for(int j=0; j<numCovs; ++j){
    arma::mat design_now = designMat1[j];
    
    design_now.each_col() %= gPi_ij;
    arma::mat grest_now= arma::sum(design_now,0);// calculat the column sum
    
    grad(span((j+1)*nk,(j+2)*nk-1 )) = trans(grest_now);
  }
  arma::vec temp2 = -2*grad;
  NumericVector gBinomLossNum = NumericVector(temp2.begin(), temp2.end());

  //---done
  
  
// Proximal gradient descent update
  List out=oneUpdateCpp(theta, current_lossval,gBinomLossNum,  intStepSize, lambda1, dat, basisMat0, nk, numCovs, shrinkScale, 
                        designMat1,theta, iter, accelrt, truncation);
  
    
  NumericVector theta_new = out["theta_l_proximal"];
  List theta_new_sep=out["theta_l_proximal_sep"];
 // stepSizeVec(iter) = out["stepSize"];
  double lossnew = out["current_lossval"];
  pi_ij_now = out["pi_ij_new"];   // the pi_ij has been updated
  
  
   // calculating stopping rules   -- criterion 1: theta between two iterations are very close
  arma::vec diff=theta_new-theta;
  double innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double tol = pow(innerdot, 0.5);
  
  // criterion 2: objective function(binomloss + penalty) are too close
  
  
  double penTerms = twoPenaltiesCpp(theta_new_sep, lambda1, numCovs, nk);
  double lossSum = lossnew + penTerms;
  double lossSumOld = lossSum;
  // 
  
  double tol2 = 100; // the second tolerance is the lossSum/observation are very close
  
  
  while(((tol > epsilon) & (iter < maxInt-1)) & (tol2 >epsilon)){
    
    lossSumOld = lossSum;
    theta = theta_new;
    current_lossval = lossnew;
    
    //--- gbinomObjectCpp ---//
    //calculate gradients
    pi_ij_arma = pi_ij_now;
    gPi_ij = meth - total % pi_ij_arma;
    
   
    
    temp=basisMat0;
    
    temp.each_col() %= gPi_ij;
    g0_now= arma::sum(temp,0); // calculat the column sum
    
    grad(span(0, nk-1)) = trans(g0_now);
    for(int j=0; j<numCovs; ++j){
      arma::mat design_now = designMat1[j];
      
       design_now.each_col() %= gPi_ij;
       arma::mat grest_now= arma::sum(design_now,0);// calculat the column sum
      
      grad(span((j+1)*nk,(j+2)*nk-1 )) = trans(grest_now);
    }
     temp2 = -2*grad;
    gBinomLossNum = NumericVector(temp2.begin(), temp2.end());
    
    //--done
//    gBinomLossNum = gbinomObjectCpp(pi_ij_now, basisMat0, dat,  nk, numCovs, designMat1);
    
    
    
    
    out=oneUpdateCpp(theta, current_lossval, gBinomLossNum,  intStepSize, lambda1, dat, basisMat0, nk, numCovs, shrinkScale, 
                          designMat1,theta, iter, accelrt, truncation);
    
    theta_new = out["theta_l_proximal"];
    theta_new_sep=out["theta_l_proximal_sep"];
    lossnew = out["current_lossval"];
    pi_ij_now = out["pi_ij_new"];   
   
    
    
    // Rcout << "step size : " << intStepSize << "\n";
    
    diff=theta_new-theta;
    innerdot = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    tol = pow(innerdot, 0.5);
    
    // Add another tolerance checking for the values of objective functions  
    
    penTerms = twoPenaltiesCpp(theta_new_sep, lambda1, numCovs, nk);
    lossSum = lossnew + penTerms;
    
    
    tol2 = lossSumOld-lossSum;
    tol2 = tol2/nsamp;
    // 
    iter = iter+1;
   // stepSizeVec(iter) = out["stepSize"];
    // arma::rowvec temp = theta_new;
    // thetaMat.row(iter) = temp;
    
    
  }
  
  
  
  List output=List::create(Named("thetaEst")=theta_new, 
                         //  Named("lossSum")=lossSum,
                          Named("neg2loglik")=lossnew,
                           Named("thetaEstSep") = theta_new_sep,
                           Named("Iter") = iter,
                          Named("gNeg2loglik")=gBinomLossNum
                          // Named("stepSizeVec") = stepSizeVec
                             );
  return(output);
  
}


