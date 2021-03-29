




#'@param tollam2  1- tollam2 is the largested value for alpha

sparseSmoothFitCV <- function(dat, n.k, stepSize=0.1, lambda = NULL, nlam = 100, lam2 = NULL, nlam2 = 10, maxInt = 500,
                              epsilon = 1E-6, printDetail = TRUE, initTheta, shrinkScale=0.5,
                              accelrt = TRUE, nfolds = 5, mc.cores,hugeCont =100000,tollam2=0.01){
  
  # shuffle the rows of dat for creating random fold for CV
  dat <- dat[sample(1:nrow(dat), nrow(dat)),]
  
  
  meth <- dat$Meth_Counts
  total <- dat$Total_Counts
  
  numCovs = ncol(dat)-4
  myp = (numCovs+1)*n.k
  lambda.min.ratio = ifelse(nrow(dat)<myp,0.01,0.0001)
  
  #---------------------------------
  # The sequence of lambda2 : ulam2
  #--------------------------------
  if (is.null(lam2)) {
    lambda_max <- 1- tollam2
    # compute lambda sequence
    #ulam2 <- seq(lambda_max, lambda_max*lambda.min.ratio, length.out = nlam2)
    
    #ulam2 <-  1-exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
    #                  length.out = nlam2))
    
    ulam2 <- seq(0, lambda_max, length.out = nlam2)
    
    #ulam2 <- c(ulam2, 0)
  } else { # user provided lambda values
    user_lambda2 = TRUE
    if (any(lam2 < 0)) stop("lambdas should be non-negative")
    ulam2 = as.double((sort(lam2)))
    nlam2 = as.integer(length(lam2))
  }
  
  
  initOut = extractMats(dat,n.k=n.k)
  #sparOmega = initOut$sparOmega
  #smoOmega1 = initOut$smoOmega1
  #designMat1 = initOut$designMat1
  #basisMat0 = initOut$basisMat0
  #-------------------------------------
  # Generate the sequence of lambda1 :
  # also need the matrix Linv
  #--------------------------------
  
  getSeqLam1HpOut = lapply(as.list(ulam2), function(x){
    getSeqLam1Hp(lambda2=x, meth, total,
                 lambda=lambda, nlam=nlam, sparOmega=initOut$sparOmega, 
                 smoOmega1=initOut$smoOmega1, designMat1=initOut$designMat1,
                 basisMat0=initOut$basisMat0, hugeCont = hugeCont)
  })
  
  lamGrid= vapply(seq(ulam2), function(i){
    getSeqLam1HpOut[[i]]$ulam
  }, FUN.VALUE =  rep(1, nlam))
  
  LinvList= lapply(seq(ulam2), function(i){
    getSeqLam1HpOut[[i]]$Linv
  })
  
  
  #--------------
  # Step 1, CV fold
  #----------------
  

  foldIndex <-caret::createFolds(meth/total, k = nfolds)
  
 # trainAll <- testAlldev <- vector("list", nfolds)
 # testAllmse <- testAlldev
  
  
  AllOut = parallel::mclapply(seq(nfolds), function(ijk){
    testID <- foldIndex[[ijk]]
    trainDat<-dat[-testID,]
    testDat <- dat[testID,]
    
    initOutTrain = extractMats(dat=trainDat,n.k=n.k)
    
    trainFit <- sparseSmoothGridRaw(meth = trainDat$Meth_Counts,
                                    total= trainDat$Total_Counts ,n.k=n.k, ulam2, lamGrid, LinvList, 
                                    theta=initTheta, stepSize=stepSize, shrinkScale=shrinkScale,
                                    basisMat0=initOutTrain$basisMat0, designMat1=initOutTrain$designMat1, 
                                    numCovs=numCovs, maxInt=maxInt, epsilon=epsilon, accelrt=accelrt,
                                    truncation=truncation, mc.cores=mc.cores)
    
    # Calculate the prediction
    
    testPred = sparseSmoothPred(trainFit=trainFit,trainDatPos=trainDat$Position,testDat=testDat,
                                 basisMat0=initOutTrain$basisMat0,basisMat1=initOutTrain$basisMat1, 
                                 n.k=n.k, numCovs=numCovs,truncation=truncation)
    
   list(testPred, trainFit$thetaOut)

  }, mc.cores=mc.cores)
  
  
  if(nlam2==1){
    testAlldev =  lapply(seq(nfolds), function(i){
      AllOut[[i]][[1]][1,]
    }) 
    
    testAllmse =  lapply(seq(nfolds), function(i){
      AllOut[[i]][[1]][2,]
    }) 
  }else{
  
  testAlldev =  lapply(seq(nfolds), function(i){
    AllOut[[i]][[1]][1,,]
  }) 
  
  testAllmse =  lapply(seq(nfolds), function(i){
    AllOut[[i]][[1]][2,,]
  }) 
}
  
 
  
  # What to export from the cv function
  # best lambda2, best lambda1
  
  #1. calculate the cross validation average/SD matrix of testPred
  
  testPredMean = Reduce("+", testAlldev) / length(testAlldev)
  
  testPredSD = apply(simplify2array(testAlldev), 1:2, sd )
  
  
  #which(testPredMean == min(testPredMean), arr.ind = TRUE)
  
  if(nlam2 >1){
  bestInd = which(testPredMean == min(testPredMean), arr.ind = TRUE)
  
  bestLambda1 = lamGrid[bestInd]
  bestLambda2 = ulam2[bestInd[2]]
  
  
  bestIndAllLam2 = apply(testPredMean, 2, which.min )
  bestLambda1vec = lamGrid[cbind(bestIndAllLam2,seq(nlam2))]
  
  }
  if(nlam2 == 1){
    bestInd = which(testPredMean == min(testPredMean), arr.ind = TRUE)
    bestLambda1 = lamGrid[bestInd,1]
    bestLambda2 = ulam2
  }
  
 
  
  #Hp1 =  (1-bestLambda2)*initOut$sparOmega + bestLambda2*initOut$smoOmega1
  #see1 = getSeqLam1HpOut[[bestInd[2]]]$Hp
  
  #initOut = extractMats(dat,n.k=n.k)
  
  if(nlam2 >1){
    
    # if a sequence of lam2 were given, report the bestFit for all lambda2
    
    
  #bestFit <- sparseSmoothPathRawOneLam1(dat[,1:3], n.k, ulam= bestLambda1, 
  #                                      lambda2=bestLambda2, Hp=getSeqLam1HpOut[[bestInd[2]]]$Hp, 
  #                                      Linv=getSeqLam1HpOut[[bestInd[2]]]$Linv, 
  #                                      theta=trainAll[[1]]$thetaOut[[bestInd[2]]][, bestInd[1]], 
  #                                      stepSize, basisMat0=initOut$basisMat0, designMat1=initOut$designMat1,
  #                                      numCovs, maxInt ,  epsilon , shrinkScale, 
  #                                      accelrt=FALSE, truncation = TRUE, basisMat1 = initOut$basisMat1)
  
  
    bestFitAll <- parallel::mclapply(seq(ulam2), function(i){
    sparseSmoothPathRawOneLam1(dat[,1:3], n.k, ulam= bestLambda1vec[i], 
                              Linv=LinvList[[i]], 
                               theta=AllOut[[1]][[2]][[i]][,bestIndAllLam2[i]], 
                               stepSize, basisMat0=initOut$basisMat0, designMat1=initOut$designMat1,
                               numCovs, maxInt ,  epsilon , shrinkScale, 
                               accelrt=FALSE, truncation = TRUE, basisMat1 = initOut$basisMat1)
  }, mc.cores=mc.cores)
   
   bestFit <- bestFitAll[[bestInd[2]]]
  
  }
  if(nlam2 == 1){
    bestFit <- sparseSmoothPathRawOneLam1(dat[,1:3], n.k, ulam= bestLambda1, 
                                          Linv=LinvList[[1]], 
                                          theta=AllOut[[1]][[2]][[1]][, bestInd[1]], 
                                          stepSize, basisMat0=initOut$basisMat0, designMat1=initOut$designMat1,
                                          numCovs, maxInt ,  epsilon , shrinkScale, 
                                          accelrt=FALSE, truncation = TRUE, basisMat1 = initOut$basisMat1)
    
  }
  
  if(nlam2==1){
  return(out = list(bestFit = bestFit,testPredMean=testPredMean,testPredSD=testPredSD,
                    ulam2=ulam2, lamGrid=lamGrid, bestLambda1=bestLambda1, bestLambda2=bestLambda2,
                    bestInd=bestInd, testAlldev=testAlldev, testAllmse=testAllmse))
  }else{
    return(out = list(bestFit = bestFit,testPredMean=testPredMean,testPredSD=testPredSD,
                      ulam2=ulam2, lamGrid=lamGrid, bestLambda1=bestLambda1, bestLambda2=bestLambda2,
                      bestInd=bestInd, bestFitAll=bestFitAll, testAlldev=testAlldev, testAllmse=testAllmse)) 
    }
}


#plotSScv <- function(cvOut){
#  testLossOut[4,,,]
#}
  
library(colortools)

plotCV <- function(cvout, hugeCont =100000, mycols = c('#a6cee3','#1f78b4','#b2df8a',
                                                       '#33a02c','#fb9a99','#e31a1c','#fdbf6f',
                                                       '#ff7f00','#cab2d6','#6a3d9a')){
  
  lamGridPlot = cvout$lamGrid
  lamGridPlot[1,] = lamGridPlot[1,]-hugeCont

  cols = sequential("steelblue", plot = FALSE)
  
  mycols <- cols[1:length(cvout$ulam2)]
  plotPred = cvout$testPredMean
  
  if(length(cvout$ulam2)>1){
  
  plot(log(lamGridPlot[,1]),  plotPred[,1], xlim = log(c(min(lamGridPlot), max(lamGridPlot))),
       ylim = c(min( plotPred), max( plotPred)), pch = 19,xlab = expression(log(lambda)), 
       ylab = "mean of CV prediction error", col = mycols[1], cex = 0.5)
  if(length(cvout$ulam2)>1){
  for(i in 2:length(cvout$ulam2)){
    points(log(lamGridPlot[,i]),  plotPred[,i], col = mycols[i], pch = 19, cex = 0.5)
  }
  }
  }
  
  if(length(cvout$ulam2)==1){
    
    plot(log(lamGridPlot[,1]),  plotPred, xlim = log(c(min(lamGridPlot), max(lamGridPlot))),
         ylim = c(min( plotPred), max( plotPred)), pch = 19,xlab = expression(log(lambda)), 
         ylab = "mean of CV prediction error", col = mycols[1], cex = 0.5)
  }
  legend('topleft', legend = paste( "alpha =", cvout$ulam2), bty = "n", cex = 0.7, fill = mycols)
  
 
}
  
#plotCV(cvout)


