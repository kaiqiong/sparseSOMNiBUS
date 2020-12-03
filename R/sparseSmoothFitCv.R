
sparseSmoothFitCV <- function(dat, n.k, stepSize=0.1,lambda1= seq(0, 5, 0.1), lambda2=seq(0, 5, 0.1), maxInt = 500,
                              epsilon = 1E-6, printDetail = TRUE, initTheta, shrinkScale=0.5,
                              accelrt = TRUE, nfolds = 10){
  dat <- dat[sample(1:nrow(dat), nrow(dat)),]
  foldIndex <-caret::createFolds(dat$Meth_Counts/dat$Total_Counts, k = nfolds)
  
  #total_index <- 1:length(foldIndex)
  #sep_index <- split(total_index, ceiling(total_index/(nfolds/ncores)))
  
  outlist <- as.list(seq(nfolds))
  for( i in seq(length(outlist))){
    outlist[[i]] <-as.list(seq(length(lambda1)*length(lambda2)))
  }
  
  gridIndex <- expand.grid(seq(length(lambda1)),seq(length(lambda2)) )
  
  testLossOut <- array(NA, dim= c(4, nfolds, length(lambda1), length(lambda2)),
                       dimnames = list( c("lam1", "lam2", "nzeros", "testLoss"), paste0("Fold", seq(nfolds)),
                                     paste0("lam1", lambda1), paste0("lam2", lambda2)))
  for ( i in seq(nfolds)){
    testID <- foldIndex[[i]]
    trainID <- setdiff(1:nrow(dat), foldIndex[[i]])
    
    trainDat<-dat[trainID,]
    testDat <- dat[testID,]
    
    initOut = extractMats(dat=trainDat,n.k=n.k)
    
     for( ll in seq(nrow(gridIndex))){
       trainFit <- fitProxGrad(theta=initTheta, stepSize=stepSize,lambda1=lambda1[gridIndex[ll,1]],
                               dat=trainDat, basisMat0= initOut$basisMat0, n.k=n.k, 
                               sparOmega= initOut$sparOmega, lambda2=lambda2[gridIndex[ll,2]], 
                               smoOmega1=initOut$smoOmega1,
                               maxInt = maxInt, epsilon = epsilon, printDetail = printDetail,shrinkScale=shrinkScale,
                               accelrt=accelrt, numCovs=initOut$numCovs, designMat1= initOut$designMat1, basisMat1=initOut$basisMat1)
       
       pred = sparseSmoothPred(fit=trainFit,dat=testDat,n.k=n.k,lambda1=lambda1[gridIndex[ll,1]],lambda2=lambda2[gridIndex[ll,2]])
       
       outlist[[i]][[ll]] = c(trainFit, pred)
       testLossOut[,i,gridIndex[ll,1],gridIndex[ll,2]]<-c(trainFit$lambda1, trainFit$lambda2, trainFit$nzeros, pred['testLoss'])
       
       pred
     }
   
       
   
  }
  
  # choose the best tune parameters
  
  # 

  return(list(outlist = outlist, testLossOut = testLossOut))
   
}


#plotSScv <- function(cvOut){
#  testLossOut[4,,,]
#}
  
  