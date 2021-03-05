#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;




double getPrec(double x);

NumericVector myseq(int &first, int &last);


List getSeparateThetaCpp(const NumericVector& theta,
                      const int& nk, 
                      const int& numCovs);


NumericVector estimatePijCpp(const List& thetaSep,
                    const arma::mat& basisMat0, 
                    const List& designMat1,
                    const int& numCovs,
                    const bool& truncation);

List binomObjectLossOnlyCpp(const NumericVector& theta,
                            const arma::mat& basisMat0, 
                            const DataFrame& dat, 
                            const int& nk,
                            const int& numCovs,
                            const List& designMat1, 
                            const bool& truncation);

double binomObjectCpp(const NumericVector& pi_ij,
                      const DataFrame& dat);

NumericVector gbinomObjectCpp(const NumericVector& pi_ij,
                              const arma::mat& basisMat0, 
                              const DataFrame& dat, 
                              const int& nk,
                              const int& numCovs,
                              const List& designMat1);
double twoPenaltiesCpp(List thetaSep,
                       double lambda1, 
                       int numCovs, 
                       int nk);
List estimatePijCpp1(const NumericVector& theta,
                     const arma::mat& basisMat0, 
                     const List& designMat1, 
                     const int& nk, 
                     const int& numCovs);