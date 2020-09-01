#' @title Prediction of a sparse-smoothness fit
#' @param fit an object from \code{sparseSmoothFit}
#' @param dat test dataset for prediction
#' @param n.k number of knots for all covariates (including intercept);
#' curretnly, we assume the same n.k for all covariates
#' @param lengthUniqueDataID number of samples in the data
#' @param numCovs number of covariates
#' @param lambda1 penalization parameter for the L2 norm
#' @param lambda2 penalization parameter for the weight between two penalities
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{var.cov.alpha} var of alpha
#' \item \code{var.alpha.0} var of alpha0
#' \item \code{var.alpha.sep} var of alpha_p, p = 1, 2, P
#' }
#' @author Kaiqiong Zhao
#' @noRd
# one theta update from the prox(theta_old - t gradient(theta_old))
sparseSmoothPred <- function(fit, dat, n.k, lambda1, lambda2){
  testMats <- extractMats(dat=dat, n.k=n.k)
  numCovs <- testMats$numCovs
  Hp <- testMats$sparOmega + lambda2*testMats$smoOmega1
  testOut <- binomObject(theta=fit$thetaEst,basisMat0=testMats$basisMat0,dat=testDat,n.k=n.k,
                         numCovs=numCovs,designMat1=testMats$designMat1)
  penTerms <- twoPenalties( getSeparateTheta(fit$thetaEst, n.k = n.k, numCovs = numCovs),
                            Hp=Hp, lambda1=lambda1, numCovs=testMats$numCovs, n.k=n.k)
  
  return(c(testLoss=testOut$neg2loglik + penTerms, binomLoss = testOut$neg2loglik,
  penTerms=penTerms))
  
}