#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


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

// [[Rcpp::export]]
List estimatePijCpp(const NumericVector& theta,
                    const arma::mat& basisMat0, 
                    const List& designMat1, 
                    const int& nk, 
                    const int& numCovs){
  
  
  List thetaSep=getSeparateThetaCpp(theta, nk, numCovs);
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
  List out=List::create(Named("pi_ij")=pi_ij, 
                        Named("theta.sep")=thetaSep);
  return(out);

      
}
// [[Rcpp::export]]
double getPrec(double x) {
  return nextafter(x, std::numeric_limits<double>::infinity()) - x;}



// [[Rcpp::export]]
List binomObjectCppValueCopy(const NumericVector theta,
                        arma::mat basisMat0, 
                         const DataFrame dat, 
                         const int nk,
                         const int numCovs,
                        List designMat1, 
                         const bool truncation){ 
  
  List estimatePijOut= estimatePijCpp(theta, basisMat0, designMat1, nk, numCovs);
  NumericVector pi_ij=estimatePijOut["pi_ij"];
  double eps=getPrec(1) * 10;
  
  if(truncation==TRUE){
    pi_ij[pi_ij > (1 - eps)] = 1 - eps;
    pi_ij[pi_ij < eps] = eps;
  }
  List thetaSep=estimatePijOut["theta.sep"];
  arma::vec meth=dat["Meth_Counts"];
  arma::vec total=dat["Total_Counts"];
  arma::vec unmeth= total-meth;
  
  arma::vec pi_ij_arma=pi_ij;
  arma::vec pi_ij_arma_c=1-pi_ij;
  arma::vec logpi=arma::log(pi_ij_arma);
  arma::vec log1mpi=arma::log(pi_ij_arma_c);
  
  arma::vec loglik_ij = meth % logpi + unmeth % log1mpi;
  
  double neg2loglik =  (-2)*arma::sum(loglik_ij);
  //calculate gradients
  
  arma::vec gPi_ij = meth - total % pi_ij_arma;
  
  int mycols=nk*(numCovs+1);
  
  arma::vec grad(mycols);
  
  basisMat0.each_col() %= gPi_ij;
  arma::mat g0_now= arma::sum(basisMat0,0);// calculat the column sum
  
  grad(span(0, nk-1)) = trans(g0_now);
 for(int j=0; j<numCovs; ++j){
   arma::mat design_now = designMat1[j];
   
   design_now.each_col() %= gPi_ij;
   arma::mat grest_now= arma::sum(design_now,0);// calculat the column sum
   
   grad(span((j+1)*nk,(j+2)*nk-1 )) = trans(grest_now);
 }
  arma::vec gNeg2loglik = -2*grad;
  List out=List::create(Named("neg2loglik")=neg2loglik, 
                         Named("theta.sep")=thetaSep,
                        Named("gNeg2loglik")=gNeg2loglik);
  
  return(out);}

// [[Rcpp::export]]
List binomObjectCpp(const NumericVector& theta,
                    const arma::mat& basisMat0, 
                    const DataFrame& dat, 
                    const int& nk,
                    const int& numCovs,
                    const List& designMat1, 
                    const bool& truncation){ 
  
  List estimatePijOut= estimatePijCpp(theta, basisMat0, designMat1, nk, numCovs);
  NumericVector pi_ij=estimatePijOut["pi_ij"];
  double eps=getPrec(1) * 10;
  
  if(truncation==TRUE){
    pi_ij[pi_ij > (1 - eps)] = 1 - eps;
    pi_ij[pi_ij < eps] = eps;
  }
  List thetaSep=estimatePijOut["theta.sep"];
  arma::vec meth=dat["Meth_Counts"];
  arma::vec total=dat["Total_Counts"];
  arma::vec unmeth= total-meth;
  
  arma::vec pi_ij_arma=pi_ij;
  arma::vec pi_ij_arma_c=1-pi_ij;
  arma::vec logpi=arma::log(pi_ij_arma);
  arma::vec log1mpi=arma::log(pi_ij_arma_c);
  
  arma::vec loglik_ij = meth % logpi + unmeth % log1mpi;
  
  double neg2loglik =  (-2)*arma::sum(loglik_ij);
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
  NumericVector gNeg2loglik =NumericVector(temp2.begin(), temp2.end());
  NumericVector exp_pi_ij = NumericVector(pi_ij.begin(), pi_ij.end());
  List out=List::create(Named("neg2loglik")=neg2loglik, 
                        Named("theta.sep")=thetaSep,
                        Named("gNeg2loglik")=gNeg2loglik,
                        Named("pi_ij")= exp_pi_ij);
  // exporting pi_ij can be deleted later
  return(out);}


// [[Rcpp::export]]
List binomObjectCppRef(const NumericVector& theta,
                    const arma::mat& basisMat0, 
                    const DataFrame& dat, 
                    const int& nk,
                    const int& numCovs,
                    const List& designMat1, 
                    const bool& truncation){ 
  
  List estimatePijOut= estimatePijCpp(theta, basisMat0, designMat1, nk, numCovs);
  NumericVector pi_ij=estimatePijOut["pi_ij"];
  double eps=getPrec(1) * 10;
  
  if(truncation==TRUE){
    pi_ij[pi_ij > (1 - eps)] = 1 - eps;
    pi_ij[pi_ij < eps] = eps;
  }
  List thetaSep=estimatePijOut["theta.sep"];
  arma::vec meth=dat["Meth_Counts"];
  arma::vec total=dat["Total_Counts"];
  arma::vec unmeth= total-meth;
  
  arma::vec pi_ij_arma=pi_ij;
  arma::vec pi_ij_arma_c=1-pi_ij;
  arma::vec logpi=arma::log(pi_ij_arma);
  arma::vec log1mpi=arma::log(pi_ij_arma_c);
  
  arma::vec loglik_ij = meth % logpi + unmeth % log1mpi;
  
  double neg2loglik =  (-2)*arma::sum(loglik_ij);
  //calculate gradients
  
  arma::vec gPi_ij = meth - total % pi_ij_arma;
  
  int myrows=basisMat0.n_rows;
  int mycols=nk*(numCovs+1);
  
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
  
  List out=List::create(Named("neg2loglik")=neg2loglik, 
                        Named("theta.sep")=thetaSep,
                        Named("gNeg2loglik")=gNeg2loglik);
  
  return(out);}



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


double lambdaMaxCpp(arma::vec& meth,
                    arma::vec& total,
                    const List& designMat1, 
                    const arma::mat& Hp, 
                    int& numCovs,
                    arma::vec& mu){
  arma::vec diffvec = meth-total % mu;
  
  arma::vec vals(numCovs);
  arma::vec Hpinv =inv(Hp);
  
  for(int i=0; i<numCovs; ++i){
    arma::mat designNow = designMat1[i];
    arma::vec bp =  2*trans(designNow) * diffvec;
    arma::rowvec temp1 =  trans(bp) * Hpinv;
    arma::vec val = temp1 * bp;
    //Rcout << "val : " << val << "\n";
    vals[i] = pow(val(0), 0.5);
  }
  
  double max_val= vals.max();
  
 // Rcout << "max : " << see << "\n";
   return(max_val);
  
}
// [[Rcpp::export]]
List designToTilda( const List& designMat1, 
                    const arma::mat& Linv,
                    const int& numCovs){ 
  List out(numCovs);
  for (int i =0; i <numCovs; ++i){
    arma::mat temp = designMat1[i];
    out[i] = temp * Linv ;
  }
  return(out);
}