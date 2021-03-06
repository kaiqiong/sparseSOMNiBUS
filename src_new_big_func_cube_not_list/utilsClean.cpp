#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;




double getPrec(double x) {
  return nextafter(x, std::numeric_limits<double>::infinity()) - x;}

NumericVector myseq(int &first, int &last) {
  NumericVector y(abs(last - first) + 1);
  if (first < last) 
    std::iota(y.begin(), y.end(), first);
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}


// [[Rcpp::export]]
List getSeparateThetaCpp1(const NumericVector& theta,
                      const int& nk, 
                      const int& numCovs){
  int size=numCovs+1;
  List thetaSep(size);
  for (int i=0; i < size; ++i){
    int lower = (i*nk);
    int upper = (i+1)*nk-1;
    thetaSep[i] = theta[myseq(lower, upper)];
  }
  return(thetaSep);
}

// estimatePijCpp directly takes the thetaSep other than theta


//'@title estimatePijCpp
//'@description return a vector of pi_ij

// [[Rcpp::export]]
NumericVector estimatePijCpp1(const List& thetaSep,
                    const arma::mat& basisMat0, 
                    const List& designMat1,
                    const int& numCovs,
                    const bool& truncation){

  int mycols=numCovs+1;
  int myrows=basisMat0.n_rows;
  
  arma::mat lpSep(myrows,mycols);
  
  arma::vec a = thetaSep[0];
  lpSep.col(0) =  basisMat0 * a; 
 // arma::mat see=lpSep;
 
  for (int i=1; i < mycols; ++i){
  
      arma::vec b = thetaSep[i];
  
      arma::mat designNow = designMat1[i-1];

     lpSep.col(i) = designNow * b;

    }
  //Rcout << "The value of sum : " << i << "\n";
  arma::vec lp_ij  = arma::sum(lpSep,1);
  //Rcout << "The value of sum : " << lp_ij(0) << "\n";
  arma::vec pi_ij = 1/(1+trunc_exp(-lp_ij));
  NumericVector pi_ij1 =NumericVector(pi_ij.begin(), pi_ij.end());
  
  double eps=getPrec(1) * 10;
  if(truncation==TRUE){
    pi_ij1[pi_ij1 > (1 - eps)] = 1 - eps;
    pi_ij1[pi_ij1 < eps] = eps;
  }
  
  return(pi_ij1);
}


//'@title binomObjectCpp
//'@description return an numeric value for neg2loglik

// [[Rcpp::export]]
double binomObjectCpp1(const NumericVector& pi_ij,
                    const DataFrame& dat){ 

  arma::vec meth=dat["Meth_Counts"];
  arma::vec total=dat["Total_Counts"];
  arma::vec unmeth= total-meth;
  
  arma::vec pi_ij_arma=pi_ij;
  arma::vec pi_ij_arma_c=1-pi_ij;
  arma::vec logpi=arma::log(pi_ij_arma);
  arma::vec log1mpi=arma::log(pi_ij_arma_c);
  
  arma::vec loglik_ij = meth % logpi + unmeth % log1mpi;
  
  double neg2loglik =  (-2)*arma::sum(loglik_ij);
  
  return(neg2loglik);
  }



//'@title gbinomObjectCpp
//'@description return an numeric vector for the gradient of the nega2loglik
// [[Rcpp::export]]
NumericVector gbinomObjectCpp1(const NumericVector& pi_ij,
                       const arma::mat& basisMat0, 
                       const DataFrame& dat, 
                       const int& nk,
                       const int& numCovs,
                       const List& designMat1){ 
  
  //calculate gradients
  arma::vec meth=dat["Meth_Counts"];
  arma::vec total=dat["Total_Counts"];
 // arma::vec unmeth= total-meth;
  
  arma::vec pi_ij_arma=pi_ij;

  
  arma::vec gPi_ij = meth - total % pi_ij_arma;
  
  int myrows=basisMat0.n_rows;
  int mycols=nk*(numCovs+1); // Revise for different nk LLLater
  
  arma::mat see(myrows, mycols);
  // this one is even slower than R function sweep
  
  for (int i=0; i < myrows ; ++i){
    arma::mat g0 =  basisMat0.row(i) * gPi_ij(i); // nk columns
    see.submat(i, 0, i, nk-1) = g0;
    for(int j=0; j<numCovs; ++j){
      arma::mat design_now = designMat1[j];
      
      arma::mat grest=design_now.row(i) * gPi_ij(i);
      
      see.submat(i,(j+1)*nk, i, (j+2)*nk-1)=grest;
    }
  }
  arma::mat grad= arma::sum(see,0);// calculat the column sum
  arma::vec gradVec=trans(grad);
  arma::vec gNeg2loglik = -2*gradVec;
  
  NumericVector gNeg2loglik1 = Rcpp::NumericVector(gNeg2loglik.begin(), gNeg2loglik.end());
  return(gNeg2loglik1);}




// [[Rcpp::export]]

double twoPenaltiesCpp(List thetaSep,
                       double lambda1, 
                       int numCovs, 
                       int nk){
  
  NumericVector squaredIndPen(numCovs);
  
  for(int i=0; i<numCovs; ++i){
    arma::rowvec thetanow = thetaSep[i+1];
    // arma::mat prod1 = thetanow * Hp;
    // Rcout << "prod1 : " << prod1 << "\n";
    // Rcout << "thetanow : " << thetanow << "\n";
    double innerdot = std::inner_product(thetanow.begin(), thetanow.end(), thetanow.begin(), 0.0);
    squaredIndPen[i] = innerdot;
  }
  double penalTerm = lambda1 * sum(sqrt(squaredIndPen ));
  
  return(penalTerm);
}



// [[Rcpp::export]]
double binomObjectCppAll (const NumericVector& theta,
   const int& nk, 
   const int& numCovs,
    const DataFrame& dat,
    const arma::mat& basisMat0, 
    const List& designMat1,
    const bool& truncation){
  List thetaSep=getSeparateThetaCpp1(theta, nk, numCovs);
  NumericVector pi_ij_now = estimatePijCpp1(thetaSep,basisMat0, designMat1, numCovs,truncation);
  
  double current_lossval = binomObjectCpp1(pi_ij_now,dat);
  
  // Calculate the gradient of -2loglik
  
 // NumericVector gBinomLossNum = gbinomObjectCpp1(pi_ij_now, basisMat0, dat,  nk, numCovs, designMat1);
  
  
  
  return(current_lossval);
}

// [[Rcpp::export]]
List getSeparateThetaCpp(const NumericVector& theta,
                          const int& nk, 
                          const int& numCovs){
  int size=numCovs+1;
  List thetaSep(size);
  for (int i=0; i < size; ++i){
    int lower = (i*nk);
    int upper = (i+1)*nk-1;
    thetaSep[i] = theta[myseq(lower, upper)];
  }
  return(thetaSep);
}

// estimatePijCpp directly takes the thetaSep other than theta


//'@title estimatePijCpp 
//'@description return a vector of pi_ij


//estimatePijCpp +  binomObjectCpp

// [[Rcpp::export]]
List binomObjectCppLossOnly(const List& thetaSep,
                             const arma::mat& basisMat0, 
                             const List& designMat1,
                             const int& numCovs,
                             const bool& truncation,
                             const DataFrame& dat){
  
  int mycols=numCovs+1;
  int myrows=basisMat0.n_rows;
  
  arma::mat lpSep(myrows,mycols);
  
  arma::vec a = thetaSep[0];
  lpSep.col(0) =  basisMat0 * a; 
  // arma::mat see=lpSep;
  
  for (int i=1; i < mycols; ++i){
    
    arma::vec b = thetaSep[i];
    
    arma::mat designNow = designMat1[i-1];
    
    lpSep.col(i) = designNow * b;
    
  }
  //Rcout << "The value of sum : " << i << "\n";
  arma::vec lp_ij  = arma::sum(lpSep,1);
  //Rcout << "The value of sum : " << lp_ij(0) << "\n";
  arma::vec pi_ij = 1/(1+trunc_exp(-lp_ij));
  NumericVector pi_ij1 =NumericVector(pi_ij.begin(), pi_ij.end());
  
  double eps=getPrec(1) * 10;
  if(truncation==TRUE){
    pi_ij1[pi_ij1 > (1 - eps)] = 1 - eps;
    pi_ij1[pi_ij1 < eps] = eps;
  }
  
  arma::vec meth=dat["Meth_Counts"];
  arma::vec total=dat["Total_Counts"];
  arma::vec unmeth= total-meth;
  
  arma::vec pi_ij_arma=pi_ij1;
  arma::vec pi_ij_arma_c=1-pi_ij1;
  arma::vec logpi=arma::log(pi_ij_arma);
  arma::vec log1mpi=arma::log(pi_ij_arma_c);
  
  arma::vec loglik_ij = meth % logpi + unmeth % log1mpi;
  
  double neg2loglik =  (-2)*arma::sum(loglik_ij);
  
  List out=List::create(Named("neg2loglik")=neg2loglik, 
                        Named("pi_ij")=pi_ij1);
  return(out);
  
}






