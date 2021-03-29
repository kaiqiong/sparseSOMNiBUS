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


double twoPenaltiesCpp(List thetaSep,
                       double lambda1, 
                       int numCovs, 
                       int nk){
  
  NumericVector squaredIndPen(numCovs);
  
  for(int i=0; i<numCovs; ++i){
    NumericVector thetanow = thetaSep[i+1];
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
//'@title estimatePijCpp 
//'@description return a vector of pi_ij


//estimatePijCpp +  binomObjectCpp


List binomObjectCppLossOnly(const List& thetaSep,
                            const arma::mat& basisMat0, 
                            const List& designMat1,
                            const int& numCovs,
                            const bool& truncation,
                            const arma::vec& meth,
                            const arma::vec& unmeth){
  
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
  
 // arma::vec meth=dat["Meth_Counts"];
 // arma::vec total=dat["Total_Counts"];
 // arma::vec unmeth= total-meth;

  arma::vec pi_ij_arma(pi_ij1.begin(), myrows, false);
  arma::vec pi_ij_arma_c=1-pi_ij_arma;
  arma::vec logpi=arma::log(pi_ij_arma);
  arma::vec log1mpi=arma::log(pi_ij_arma_c);
  
  arma::vec loglik_ij = meth % logpi + unmeth % log1mpi;
  
  double neg2loglik =  (-2)*arma::sum(loglik_ij);
  
  List out=List::create(Named("neg2loglik")=neg2loglik, 
                        Named("pi_ij")=pi_ij_arma);
  return(out);
  
}



// [[Rcpp::export]]
List binomObjectCppLossOnlyVec(const NumericVector& theta,
                            const arma::mat& basisMat0, 
                            const int& numCovs,
                            const List& designMat1, 
                            const bool& truncation,
                            const arma::vec& meth,
                            const arma::vec& unmeth,
                            const arma::vec& total,
                            const int& nk){ 
  List thetaSep = getSeparateThetaCpp(theta, nk, numCovs);
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
  
  // arma::vec meth=dat["Meth_Counts"];
  // arma::vec total=dat["Total_Counts"];
  // arma::vec unmeth= total-meth;
  
  arma::vec pi_ij_arma(pi_ij1.begin(), myrows, false);
  arma::vec pi_ij_arma_c=1-pi_ij_arma;
  arma::vec logpi=arma::log(pi_ij_arma);
  arma::vec log1mpi=arma::log(pi_ij_arma_c);
  
  arma::vec loglik_ij = meth % logpi + unmeth % log1mpi;
  
  double neg2loglik =  (-2)*arma::sum(loglik_ij);
  
  

  arma::vec temp=meth/total;
  arma::vec temp2 = temp-pi_ij_arma;
  double res1 = std::inner_product(temp2.begin(), temp2.end(), temp2.begin(), 0.0);
  // double res=sqrt(dot(temp1, u_p));
  double res=sqrt(res1);

  
  List out=List::create(Named("neg2loglik")=neg2loglik, 
                        Named("sqrt_of_sum_of_dif")=res);
  return(out);
  
  }


//'@title gbinomObjectCpp
//'@description return an numeric vector for the gradient of the nega2loglik
// [[Rcpp::export]]
NumericVector gbinomObjectCpp(const arma::vec& pi_ij_arma,
                               const arma::mat& basisMat0, 
                               const arma::vec& meth,
                               const arma::vec& total,
                               const int& nk,
                               const int& numCovs,
                               const List& designMat1){ 
  
  
  //calculate gradients
  
  int mycols=nk*(numCovs+1);
  arma::vec gPi_ij = meth - total % pi_ij_arma;
  
  
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

  
  return(gBinomLossNum);}


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


//'@title thetaUpdateCpp
//'@description One proximal gradient descent update given the stepSize, current theta, and current gradient

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

List oneUpdateCpp(const NumericVector& theta,
                  const double& current_lossval,
                  const NumericVector& gBinomLossNum,
                  double stepSize,    // Not adding the pointer for stepSize variable is important, other wise, stepSize will be changed after running the line searching
                  const double& lambda1,
                  const arma::vec& meth,
                  const arma::vec& unmeth,
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
  
  //  NumericVector pi_ij_now = estimatePijCpp(theta_l_p_sep,basisMat0, designMat1, numCovs,truncation);
  //   double newloss = binomObjectCpp(pi_ij_now,dat);
  
  List out = binomObjectCppLossOnly(theta_l_p_sep,basisMat0, designMat1, numCovs,truncation, meth, unmeth);
  double newloss = out["neg2loglik"];
  
  
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
    
    
    
    //  pi_ij_now = estimatePijCpp(theta_l_p_sep,basisMat0, designMat1, numCovs,truncation);
    //  newloss = binomObjectCpp(pi_ij_now,dat);
    
    
    out = binomObjectCppLossOnly(theta_l_p_sep,basisMat0, designMat1, numCovs,truncation, meth, unmeth);
    newloss = out["neg2loglik"];
    //--------------------//
    
    double innerdot1 = std::inner_product(gBinomLossNum.begin(), gBinomLossNum.end(), Gttheta.begin(), 0.0);
    double innerdot2 = std::inner_product(Gttheta.begin(), Gttheta.end(), Gttheta.begin(), 0.0);
    linearSupport =  current_lossval - stepSize * innerdot1 + stepSize/2*innerdot2;
    
    
    shrinkCondi = newloss >linearSupport;
  }
  
  
  
  NumericVector newtheta=new_out["theta_l_proximal"];
  NumericVector pi_ij_now = out["pi_ij"];
  List output=List::create(Named("theta_l_proximal")=newtheta, 
                           Named("stepSize")=stepSize,
                           Named("theta_l_proximal_sep")= theta_l_p_sep,
                           Named("pi_ij_new") = pi_ij_now,
                           Named("current_lossval")=newloss);
  return(output);
  
}




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
List fitProxGradCppClean1(NumericVector& theta, 
                         double& intStepSize,
                         const double& lambda1,
                         const arma::vec& meth,
                         const arma::vec& unmeth,
                         const arma::vec& total,
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
  int nsamp = meth.n_rows;
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
  
  
  List binomOut = binomObjectCppLossOnly(thetaSep, basisMat0,  designMat1, numCovs, truncation, meth, unmeth);
  
  double current_lossval =  binomOut["neg2loglik"];
  arma::vec pi_ij_arma = binomOut["pi_ij"];
  
  
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
  List out=oneUpdateCpp(theta, current_lossval,gBinomLossNum,  intStepSize, lambda1, meth, unmeth, basisMat0, nk, numCovs, shrinkScale, 
                        designMat1,theta, iter, accelrt, truncation);
  
  
  NumericVector theta_new = out["theta_l_proximal"];
  List theta_new_sep=out["theta_l_proximal_sep"];
  // stepSizeVec(iter) = out["stepSize"];
  double lossnew = out["current_lossval"];
  
  arma::vec seenow = out["pi_ij_new"];   // the pi_ij has been updated
  pi_ij_arma=seenow;
  
 // arma::vec pi_ij_arma(seenow.begin(), nsamp , false);
 // pi_ij_arma = out["pi_ij_new"];   // the pi_ij has been updated
  
  
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
  //  pi_ij_arma = pi_ij_now;
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
    
    
    
    
    out=oneUpdateCpp(theta, current_lossval, gBinomLossNum,  intStepSize, lambda1, meth, unmeth, basisMat0, nk, numCovs, shrinkScale, 
                     designMat1,theta, iter, accelrt, truncation);
    
    theta_new = out["theta_l_proximal"];
    theta_new_sep=out["theta_l_proximal_sep"];
    lossnew = out["current_lossval"];

    arma::vec seenow1 = out["pi_ij_new"];   
    pi_ij_arma=seenow1;
    
    
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


// I removed the PathRaw and GridRaw, which are now inside the folder 'src_new_big_func'
