#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


NumericVector myseq(int &first, int &last);

List getSeparateThetaCpp(const NumericVector& theta,
                      const int& nk, 
                      const int& numCovs);

List estimatePijCpp(const NumericVector& theta,
                    const arma::mat& basisMat0, 
                    const List& designMat1, 
                    const int& nk, 
                    const int& numCovs);

double getPrec(double x);



List binomObjectCpp(const NumericVector& theta,
                    const arma::mat& basisMat0, 
                    const DataFrame& dat, 
                    const int& nk,
                    const int& numCovs,
                    const List& designMat1, 
                    const bool& truncation);

List binomObjectLossOnlyCpp(const NumericVector& theta,
                            const arma::mat& basisMat0, 
                            const DataFrame& dat, 
                            const int& nk,
                            const int& numCovs,
                            const List& designMat1, 
                            const bool& truncation);

double twoPenaltiesCpp(List thetaSep,
                       double lambda1, 
                       int numCovs, 
                       int nk);
double lambdaMaxCpp(const DataFrame& dat, 
                    const List& designMat1, 
                    const arma::mat& Hp, 
                    int& numCovs);