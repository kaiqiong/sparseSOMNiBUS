#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "utilsClean.h"



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



//'@title thetaUpdateCpp
//'@description One proximal gradient descent update given the stepSize, current theta, and current gradient



// [[Rcpp::export]]
List thetaUpdateCpp(const double& stepSize,
                    const NumericVector& theta, 
                    const NumericVector& gBinomLoss, 
                    const int& nk,
                    const int& numCovs, 
                    const double& lambda1, 
                    const NumericVector& theta_m,
                    const int& iter,
                    const bool& accelrt){
  // delete the acceleration part
  //if(accelrt == TRUE){
  // NumericVector thetaInterm = theta + ((iter - 2)/(iter +1)) *(theta-theta_m);
  //  List  binom_outInterm = binomObjectCpp(thetaInterm,basisMat0, dat, nk,numCovs, designMat1, truncation);
  //  NumericVector theta_l =  thetaInterm - stepSize * binom_outInterm["gNeg2loglik"];
  //}else{
  NumericVector theta_l =  theta - stepSize * gBinomLoss;
  //}
  List thetaL_sep = getSeparateThetaCpp(theta_l, nk, numCovs);
  List theta_l_p_sep(numCovs+1);
  
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
    theta_l_p_sep[i] = temp;
    theta_l_proximal[myseq(lower, upper)] = temp ;
  } 
  NumericVector  G_t_theta = (theta - theta_l_proximal)/stepSize;
  

  
  List out=List::create(Named("G_t_theta")= G_t_theta, 
                        Named("theta_l_proximal")=theta_l_proximal,
                        Named("theta_l_p_sep")=theta_l_p_sep); 
  
  return(out);
  
}


//'@title oneUpdateCpp
//'@description One proximal gradient descent update given the stepSize


// [[Rcpp::export]]
List oneUpdateCpp(const NumericVector& theta,
                  const double& current_lossval,
                  const NumericVector& gBinomLossNum,
               double stepSize,    // Not adding the pointer for stepSize variable is important, other wise, stepSize will be changed after running the line searching
               const double& lambda1, 
               const DataFrame& dat,
               const arma::mat& basisMat0,
               const int& nk,
               const int& numCovs,
               const double& shrinkScale,
               const List& designMat1, 
               const NumericVector& theta_m,  // useless argument put here if we want to implement accelerationg later on
               const int& iter, 
               const bool& accelrt,
               const bool& truncation){

List new_out = thetaUpdateCpp(stepSize,theta, gBinomLossNum,
                          nk, numCovs, lambda1, theta_m, iter, accelrt);


//if(accelrt ==TRUE){
//  arma::vec tempDiff = new_out["theta_l_proximal"] - new_out["thetaInterm"];
//  double linearSupport =  new_out["binom_outInterm"]["neg2loglik"] + dot(new_out["binom_outInterm"]["neg2loglik"], tempDiff) +
//    dot(tempDiff, tempDiff)/(2*stepSize);
//}else{

  NumericVector Gttheta = new_out["G_t_theta"];
  List theta_l_p_sep=new_out("theta_l_p_sep");
  
  // calculate the lossval under the new theta_l_proximal --- this is the left hand side of checking condition
  
  // NumericVector pi_ij_now = estimatePijCpp(theta_l_p_sep,basisMat0, designMat1, numCovs,truncation);
  //  double newloss = binomObjectCpp(pi_ij_now,dat);
  
  
  
  // -- get rid of calling functions here
  
  //---estimatePijCpp
  int myrows=basisMat0.n_rows;
  int size=numCovs+1;
  arma::mat lpSep(myrows,size);
  
  arma::vec a = theta_l_p_sep[0];
  lpSep.col(0) =  basisMat0 * a; 
  
  for (int i=1; i < size; ++i){
    
    arma::vec b = theta_l_p_sep[i];
    
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
  
  //---binomObjectCpp 
  arma::vec meth=dat["Meth_Counts"];
  arma::vec total=dat["Total_Counts"];
  arma::vec unmeth= total-meth;
  
  arma::vec pi_ij_arma=pi_ij_now;
  arma::vec pi_ij_arma_c=1-pi_ij_now;
  arma::vec logpi=arma::log(pi_ij_arma);
  arma::vec log1mpi=arma::log(pi_ij_arma_c);
  
  arma::vec loglik_ij = meth % logpi + unmeth % log1mpi;
  
  double newloss =  (-2)*arma::sum(loglik_ij);
  
  //--done
 
  
  double innerdot1 = std::inner_product(gBinomLossNum.begin(), gBinomLossNum.end(), Gttheta.begin(), 0.0);
  double innerdot2 = std::inner_product(Gttheta.begin(), Gttheta.end(), Gttheta.begin(), 0.0);
  double linearSupport =  current_lossval - stepSize * innerdot1 + stepSize/2*innerdot2;
//}

  bool   shrinkCondi = newloss >linearSupport;
    while(shrinkCondi == TRUE){
      stepSize=stepSize*shrinkScale;
      new_out = thetaUpdateCpp(stepSize,theta, gBinomLossNum,
                               nk, numCovs, lambda1, theta_m, iter, accelrt);
      Gttheta = new_out["G_t_theta"];
      
      theta_l_p_sep=new_out("theta_l_p_sep");
      
      
      
      //----- Expand this two lines of code
      
      
     // pi_ij_now = estimatePijCpp(theta_l_p_sep,basisMat0, designMat1, numCovs,truncation);
    //  newloss = binomObjectCpp(pi_ij_now,dat);
      
      // -- ----------------//
      
      //---estimatePijCpp
    
      
      arma::vec aa = theta_l_p_sep[0];
      lpSep.col(0) =  basisMat0 * aa; 
      
      for (int i=1; i < size; ++i){
        
        arma::vec bb = theta_l_p_sep[i];
        
        arma::mat designNow = designMat1[i-1];
        
        lpSep.col(i) = designNow * bb;
        
      }
      
       lp_ij  = arma::sum(lpSep,1);
      
       pi_ij = 1/(1+trunc_exp(-lp_ij));
      
        pi_ij_now =NumericVector(pi_ij.begin(), pi_ij.end());
      
      if(truncation==TRUE){
        pi_ij_now[ pi_ij_now > (1 - eps)] = 1 - eps;
        pi_ij_now[ pi_ij_now < eps] = eps;
      }
      
      //---done
      
      //---binomObjectCpp 
      
      pi_ij_arma=pi_ij_now;
      pi_ij_arma_c=1-pi_ij_now;
      logpi=arma::log(pi_ij_arma);
      log1mpi=arma::log(pi_ij_arma_c);
      
      loglik_ij = meth % logpi + unmeth % log1mpi;
      
      double newloss =  (-2)*arma::sum(loglik_ij);
      
      //--done
      
      
      //--------------------//
      
      double innerdot1 = std::inner_product(gBinomLossNum.begin(), gBinomLossNum.end(), Gttheta.begin(), 0.0);
      double innerdot2 = std::inner_product(Gttheta.begin(), Gttheta.end(), Gttheta.begin(), 0.0);
      linearSupport =  current_lossval - stepSize * innerdot1 + stepSize/2*innerdot2;
     
     
      shrinkCondi = newloss >linearSupport;
    }
    
    
    
    NumericVector newtheta=new_out["theta_l_proximal"];

    List out=List::create(Named("theta_l_proximal")=newtheta, 
                          Named("stepSize")=stepSize,
                          Named("theta_l_proximal_sep")= theta_l_p_sep,
                          Named("pi_ij_new") = pi_ij_now,
                          Named("current_lossval")=newloss);
    return(out);

}


