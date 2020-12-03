#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "utils.h"


template <typename T>
NumericVector arma2vec(const T& x){
  return Rcpp::NumericVector(x.begin(), x.end());
}

// [[Rcpp::export]]
NumericVector proximalOperatorCpp(const double& t, 
                                  const double& lambda1, 
                                  const arma::colvec& u_p,
                                  const int& nk){
  
  //NumericVector see1=NumericVector(temp1.begin(), temp2.end());
 // NumericVector see2=NumericVector(u_p.begin(), u_p.end());
 double res1 = std::inner_product(u_p.begin(), u_p.end(), u_p.begin(), 0.0);
 // double res=sqrt(dot(temp1, u_p));
 double res=sqrt(res1);
  double threshold =1-(t*lambda1)/res;
  if(threshold >=0){
    arma::colvec out =threshold *u_p;
    return(Rcpp::NumericVector(out.begin(), out.end()));
  }else{
    arma::colvec out(nk);
    out.fill(0);
    return(Rcpp::NumericVector(out.begin(), out.end()));
  }
}

//[[Rcpp::export]]

double dot_std(NumericVector& x,
               NumericVector& y){
  double res1 = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
  return(res1);
}

//[[Rcpp::export]]

double dot_arma(arma::vec& x,
                arma::vec& y){
  double res1 = dot(x,y);
  return(res1);
}
// [[Rcpp::export]]
NumericVector proximalOperatorCppOld(const double& t, 
                                  const double& lambda1, 
                                  const arma::mat& Hp, 
                                  const arma::colvec& u_p,
                                  const int& nk){
  arma::mat temp1=  Hp*u_p;
  double res=sqrt(dot(temp1, u_p));
  double threshold =1-(t*lambda1)/res;
  if(threshold >=0){
    arma::colvec out =threshold *u_p;
    return(Rcpp::NumericVector(out.begin(), out.end()));
  }else{
    arma::colvec out(nk);
    out.fill(0);
    return(Rcpp::NumericVector(out.begin(), out.end()));
  }
}


// [[Rcpp::export]]
List thetaUpdateCpp(const double& stepSize,
                    const NumericVector& theta, 
                    const NumericVector& gBinomLoss, 
                    const int& nk,
                    const int& numCovs, 
                    const double& lambda1, 
                    const arma::mat& Hp, 
                    const arma::mat& basisMat0, 
                    const DataFrame& dat, 
                    const List& designMat1,
                    const NumericVector& theta_m,
                    const int& iter,
                    const bool& accelrt,
                    const bool& truncation){
  // delete the acceleration part
  //if(accelrt == TRUE){
  // NumericVector thetaInterm = theta + ((iter - 2)/(iter +1)) *(theta-theta_m);
  //  List  binom_outInterm = binomObjectCpp(thetaInterm,basisMat0, dat, nk,numCovs, designMat1, truncation);
  //  NumericVector theta_l =  thetaInterm - stepSize * binom_outInterm["gNeg2loglik"];
  //}else{
  NumericVector theta_l =  theta - stepSize * gBinomLoss;
  //}
  List thetaL_sep = getSeparateThetaCpp(theta_l, nk, numCovs);
  //List theta_l_p_sep(numCovs+1);
  
  NumericVector theta_l_proximal(theta.length());
  
  for (int i=0; i < numCovs+1; ++i){
    NumericVector temp;
    if( i ==0){
      temp= thetaL_sep[i];
    }else{
      temp=proximalOperatorCpp(stepSize, lambda1, thetaL_sep[i], nk);
    }
    
    int lower = (i*nk);
    int upper = (i+1)*nk-1;
  //  theta_l_p_sep[i] = temp;
    theta_l_proximal[myseq(lower, upper)] = temp ;
  } 
  NumericVector  G_t_theta = (theta - theta_l_proximal)/stepSize;
  
  List binom_out_new = binomObjectCpp(theta_l_proximal,basisMat0, dat, nk, numCovs, designMat1, truncation);
  
  List out=List::create(Named("binom_out_new")=binom_out_new,
                        Named("G_t_theta")= G_t_theta, 
                        Named("theta_l_proximal")=theta_l_proximal); 
  
  return(out);
  
}


// [[Rcpp::export]]

List oneUpdateCpp(const NumericVector& theta,
               double stepSize,    // Not adding the pointer for stepSize variable is important, other wise, stepSize will be changed after running the line searching
               const double& lambda1, 
               const DataFrame& dat,
               const arma::mat& basisMat0,
               const int& nk,
               const arma::mat& Hp, 
               const int& numCovs,
               const double& shrinkScale,
               const List& designMat1, 
               const NumericVector& theta_m, 
               const int& iter, 
               const bool& accelrt,
               const bool& truncation){
List binom_out = binomObjectCpp(theta,basisMat0, dat, nk,numCovs, designMat1, truncation);
arma::vec  gBinomLoss = binom_out["gNeg2loglik"];
NumericVector gBinomLossNum = Rcpp::NumericVector(gBinomLoss.begin(), gBinomLoss.end());
List new_out = thetaUpdateCpp(stepSize,theta, gBinomLossNum,
                          nk, numCovs, lambda1, Hp,
                          basisMat0, dat,designMat1,theta_m,
                       iter, accelrt,truncation);
//if(accelrt ==TRUE){
//  arma::vec tempDiff = new_out["theta_l_proximal"] - new_out["thetaInterm"];
//  double linearSupport =  new_out["binom_outInterm"]["neg2loglik"] + dot(new_out["binom_outInterm"]["neg2loglik"], tempDiff) +
//    dot(tempDiff, tempDiff)/(2*stepSize);
//}else{
  NumericVector Gttheta = new_out["G_t_theta"];
  double lossval = binom_out["neg2loglik"];
  double innerdot1 = std::inner_product(gBinomLoss.begin(), gBinomLoss.end(), Gttheta.begin(), 0.0);
  double innerdot2 = std::inner_product(Gttheta.begin(), Gttheta.end(), Gttheta.begin(), 0.0);
  double linearSupport =  lossval - stepSize * innerdot1 + stepSize/2*innerdot2;
//}
  List newLossOut = new_out["binom_out_new"];
  double newloss = newLossOut["neg2loglik"];
  bool   shrinkCondi = newloss >linearSupport;
    while(shrinkCondi == TRUE){
      stepSize=stepSize*shrinkScale;
      new_out = thetaUpdateCpp(stepSize,theta, gBinomLossNum,
                               nk, numCovs, lambda1, Hp,
                               basisMat0, dat,designMat1,theta_m,
                               iter, accelrt,truncation);
      Gttheta = new_out["G_t_theta"];
      double innerdot1 = std::inner_product(gBinomLoss.begin(), gBinomLoss.end(), Gttheta.begin(), 0.0);
      double innerdot2 = std::inner_product(Gttheta.begin(), Gttheta.end(), Gttheta.begin(), 0.0);
      linearSupport =  lossval - stepSize * innerdot1 + stepSize/2*innerdot2;
      newLossOut = new_out["binom_out_new"];
      newloss = newLossOut["neg2loglik"];
      shrinkCondi = newloss >linearSupport;
    }
    
    NumericVector newtheta=new_out["theta_l_proximal"];
    List newthetaSep=newLossOut["theta.sep"];
    List out=List::create(Named("G_t_theta")=Gttheta, 
                          Named("theta_l_proximal")=newtheta, 
                          Named("binom_out_new")=newLossOut,
                          Named("stepSize")=stepSize,
                          Named("theta_l_proximal_sep")= newthetaSep);
    return(out);

}


