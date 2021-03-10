

# Given a fixed of grid of lam and lam2

sparseSmoothGridRaw <- function(meth, total, n.k, ulam2, lamGrid, LinvList,theta, stepSize, shrinkScale,
                             basisMat0,  designMat1, numCovs,
                             maxInt = 10^5,  epsilon = 1E-20, accelrt=FALSE, truncation = TRUE,
                             mc.cores = 1){
 
  nlam2 = length(ulam2)
  unmeth = total - meth
  neg2loglikSat <- getSatNeg2Loglik(meth, total)
  
  AllOut = lapply(seq(nlam2), function(i){
    sparseSmoothPathRaw(meth, unmeth,total, n.k, ulam=lamGrid[,i],
                        Linv=LinvList[[i]], theta, stepSize, basisMat0, 
                        designMat1,  numCovs, maxInt,  epsilon, shrinkScale,
                        accelrt, truncation, neg2loglikSat)
  })
  
  
  thetaOut<-  lapply(seq(nlam2), function(i){ Matrix::Matrix(AllOut[[i]]$thetaMat,sparse = TRUE)})
  thetaOutOri <- lapply(seq(nlam2), function(i){Matrix::Matrix(AllOut[[i]]$thetaMatOriginal,sparse = TRUE)})
  zeroCovsBool <- lapply(seq(nlam2), function(i){AllOut[[i]]$zeroCovsBool})
  IterNum <- vapply(seq(nlam2), function(i){AllOut[[i]]$Iter}, FUN.VALUE = rep(0, nrow(lamGrid)))
  
  return(out = list(thetaOut=thetaOut, 
                    thetaOutOri=thetaOutOri, 
                    zeroCovsBool=zeroCovsBool,
                    IterNum=IterNum))
}

# Given the sequence of lambda1 -- ulam
# valuesl of lambda2, Hp, Linv
# 

# Raw function takes a vector of ulam and a value of lambda2 -- better for cv


#'@param y meth_counts
#'@param x total_counts
#'@param designMat1  design matrix in the raw scale
#'@param basisMat0 design matrix in the raw scale
#'@param ulam a sequence of lambda
#'@param Linv corresponding Linv for the given lambda2
sparseSmoothPathRaw <- function(meth, unmeth, total,n.k, ulam,  Linv, theta, stepSize, basisMat0, 
                             designMat1,  numCovs, maxInt = 10^5,  epsilon = 1E-20, shrinkScale,
                             accelrt=FALSE, truncation = TRUE,neg2loglikSat){
  
  basisMat0_tilda <- basisMat0 %*% Linv 
  designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})
  
  # calculate the saturated
  
  myp = (numCovs+1)*n.k

  nlam = as.integer(length(ulam))
 
  
  thetaMat <- matrix(NA, nrow = myp, ncol = nlam )
  zeroCovsBool <- matrix(NA, nrow = numCovs, ncol = nlam)

 
  
  Iter <- rep(NA, nlam)

  fit1 <- fitProxGradCppClean1(theta, intStepSize = stepSize, lambda1 = ulam[1], meth, unmeth, total, basisMat0_tilda, n.k,
                         maxInt, epsilon, shrinkScale,
                         accelrt, numCovs, designMat1_tilda, truncation)
 

  thetaMat[,1] <- fit1$thetaEst
 
  
  zeroCovsBool[,1] <-unlist(lapply(fit1$thetaEstSep[-1], function(x){all(x==0)}))
  
  Iter[1] <- fit1$Iter
  
  for( i in 2:length(ulam)){
    
    fit1 <- fitProxGradCppClean1(fit1$thetaEst, intStepSize = stepSize, lambda1 = ulam[i], meth,unmeth, total, basisMat0_tilda, n.k,
                           maxInt, epsilon, shrinkScale,
                           accelrt, numCovs, designMat1_tilda, truncation)
    thetaMat[,i] <- fit1$thetaEst
   
    # checkall[,i] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[i], Hp, L,Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
    
    zeroCovsBool[,i] <-unlist(lapply(fit1$thetaEstSep[-1], function(x){all(x==0)}))
    Iter[i] <- fit1$Iter
    if(fit1$neg2loglik < neg2loglikSat) break
    
  }
  
  thetaMatOriginal <- 
    vapply(seq(ulam), function(i){
      unlist(lapply(getSeparateThetaCpp(thetaMat[,i], n.k, numCovs), function(x){Linv%*%x}))
    }, FUN.VALUE = rep(1, myp))
  
  return(out = list(thetaMat=thetaMat,  
                    thetaMatOriginal=thetaMatOriginal,
                    zeroCovsBool=zeroCovsBool, Iter = Iter))
  
}



sparseSmoothPathRawOneLam1 <- function(dat, n.k, ulam, Linv, theta, stepSize, basisMat0, 
                                designMat1,  numCovs, maxInt = 10^5,  epsilon = 1E-20, shrinkScale,
                                accelrt=FALSE, truncation = TRUE, basisMat1){
  
  
  basisMat0_tilda <- basisMat0 %*% Linv 
  designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})
  
  # calculate the saturated
  meth = dat$Meth_Counts
  total = dat$Total_Counts
  unmeth = total - meth
  

  fit1 <- fitProxGradCppClean1(theta, intStepSize = stepSize, lambda1 = ulam[1], meth,unmeth, total, basisMat0_tilda, n.k,
                       maxInt, epsilon, shrinkScale,
                       accelrt, numCovs, designMat1_tilda, truncation)
  

  thetaMat <- fit1$thetaEst
  
  zeroCovsBool<-unlist(lapply(fit1$thetaEstSep[-1], function(x){all(x==0)}))
  
  Iter <- fit1$Iter
  #fit1$thetaEstSep
  #checkall[,1] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[1], Hp, L, Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )

  thetaMatOriginal <- unlist(lapply(getSeparateThetaCpp(thetaMat, n.k, numCovs), function(x){Linv%*%x}))

  uniqPos = unique(dat$Position)
  
  uni_rows <- match(uniqPos, dat$Position)
  
  #penalBetas <- basisMat1[uni_rows,]%*% thetaMatOriginal
  
  thetaMatOriginalSep <- getSeparateThetaCpp(thetaMatOriginal, n.k, numCovs)
  penalBetas <-  lapply(thetaMatOriginalSep[-1][!zeroCovsBool], function(x){basisMat1[uni_rows,]%*% x})
  
  return(out = list(thetaMat=Matrix::Matrix(thetaMat, sparse=TRUE), ulam= ulam, 
                    thetaMatOriginal=Matrix::Matrix(thetaMatOriginal, sparse=TRUE),
                    zeroCovsBool=zeroCovsBool, Iter = Iter,
                    penalBetas = penalBetas, uniqPos = uniqPos))
  
}







