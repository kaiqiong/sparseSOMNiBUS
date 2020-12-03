#' @title Prediction of a sparse-smoothness fit
#' @param fit an object from \code{sparseSmoothFit} or \code{fitProxGrad}
#' @param dat test dataset for prediction
#' @param n.k number of knots for all covariates (including intercept);
#' curretnly, we assume the same n.k for all covariates
#' @param lengthUniqueDataID number of samples in the data
#' @param numCovs number of covariates
#' @param lambda1 penalization parameter for the L2 norm
#' @param lambda2 penalization parameter for the weight between two penalities
#' @return This function return a vector including objects:
#' \itemize{
#' \item \code{testLoss} var of alpha
#' \item \code{binomLoss} var of alpha0
#' \item \code{penTerms} var of alpha_p, p = 1, 2, P
#' }
#' @author Kaiqiong Zhao
#' @noRd
# one theta update from the prox(theta_old - t gradient(theta_old))
sparseSmoothPredOld <- function(fit, dat, n.k, lambda1, lambda2){
  
  rowsID <- match(dat$Position, fit$uniPos)
  if (any(is.na(rowsID))){message(paste0("Some positions in the test dataset are not present in the train fit; predictions at those postions are not available"))}
  basisMat0 <- fit$basisMat0[rowsID,]
  basisMat1 <- fit$basisMat1[rowsID,]
  designMat1 <- extractDesignMat1(numCovs=fit$numCovs, basisMat1, dat)
  #nrow(basisMat1)== nrow(dat) 
  #testMats <- extractMats(dat=dat, n.k=n.k)
  #numCovs <- testMats$numCovs
  Hp <- fit$sparOmega + lambda2*fit$smoOmega1
  testOut <- binomObject(theta=fit$thetaEst,basisMat0=basisMat0,dat = dat,n.k=n.k,
                         numCovs=fit$numCovs,designMat1=designMat1)
  penTerms <- twoPenalties( getSeparateTheta(fit$thetaEst, n.k = n.k, numCovs = fit$numCovs),
                            Hp=Hp, lambda1=lambda1, numCovs=fit$numCovs, n.k=n.k)
  
  return(c(testLoss=testOut$neg2loglik + penTerms, binomLoss = testOut$neg2loglik,
  penTerms=penTerms))
  
}



# Calculate the loss for a grid of sparse-smoothness fit

#' @title Prediction of a sparse-smoothness fit
#' @param trainFit an object from \code{sparseSmoothGrid}
#' @param trainDatPos the position vector in the trainDat \code{out}
#' @param testDat the test data set for the calculation of Loss, requires all the columns, meth, total, pos, all covariates
#' @param n.k number of knots for all covariates (including intercept);
#' @param basisMat0 basisMat0 for the trainDat
#' @param basisMat1
#' curretnly, we assume the same n.k for all covariates
#' @return This function return a vector including objects:
sparseSmoothPred <- function(trainFit, trainDatPos, testDat, basisMat0, basisMat1, n.k, numCovs, truncation){
  
  rowsID <- match(testDat$Position, unique(trainDatPos))
  if (any(is.na(rowsID))){message(paste0("Some positions in the test dataset are not present in the train fit;
                                         predictions at those postions are not available"))}
  # calculate the design matrix for the test dataset
  basisMat0 <- basisMat0[rowsID,]
  basisMat1 <- basisMat1[rowsID,]
  designMat1 <- extractDesignMat1(numCovs, basisMat1, testDat)
  
  # calculate criterion
  # 1, loss function value / number of observations
  # 2, mean prediction errors
  
  
  nlam2 = length(trainFit$thetaOut)
  nlam = lapply(trainFit$thetaOut, ncol)
  
  lossvals <- 
  vapply(seq(nlam2), function(i){
    vapply(seq(nlam[[i]]), function(j){
      
      testOut <-binomObjectCpp(theta=trainFit$thetaOutOri[[i]][,j], basisMat0=basisMat0,dat=testDat[,1:2],nk=n.k,
                               numCovs=,designMat1=designMat1, truncation=truncation 
                               )
      testOut$neg2loglik/nrow(testDat)
    }, FUN.VALUE = 1)
  }, FUN.VALUE = rep(1, nlam[[1]]))
  
 return(lossvals)
  
}


