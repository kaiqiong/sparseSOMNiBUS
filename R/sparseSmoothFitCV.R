





#'@param ulam2 a sequence of lambda2
#'
getSeqLam1Hp <- function(lambda2,dat, lambda=NULL, nlam, sparOmega, smoOmega1, designMat1, basisMat0, hugeCont =100000 ){
  Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0
  
  #--- Matrix decomposition for Hp and calculate transformed design matrix
  
  L  = chol(Hp)
  Hinv = chol2inv(L)
  Linv = solve(L)
  #Linv = MASS::ginv(L)
  
  start_fit <-  getStart(y=dat$Meth_Counts, x=dat$Total_Counts, 
                         designMat1=designMat1 , 
                         Hp=Hp,Hp_inv=Hinv, numCovs=numCovs, basisMat0=basisMat0)
  
  myp = (numCovs+1)*n.k
  lambda.min.ratio = ifelse(nrow(dat)<myp,0.01,0.0001)
  if (is.null(lambda)) {
    if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
    
    # compute lambda max: to add code here
    lambda_max <- start_fit$lambda_max
    
    # compute lambda sequence
    ulam <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                    length.out = nlam))
  } else { # user provided lambda values
    user_lambda = TRUE
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam = as.double(rev(sort(lambda)))
    nlam = as.integer(length(lambda))
  }
  
  ulam[1] = ulam[1]+hugeCont
  return(out = list(Hp=Hp, Linv=Linv, start_fit=start_fit, ulam = ulam, nlam = nlam,
                    L=L, Hinv =Hinv ))
  
}




sparseSmoothFitCV <- function(dat, n.k, stepSize=0.1, lambda = NULL, nlam = 100, lam2 = NULL, nlam2 = 10, maxInt = 500,
                              epsilon = 1E-6, printDetail = TRUE, initTheta, shrinkScale=0.5,
                              accelrt = TRUE, nfolds = 5, mc.cores){
  
  numCovs = ncol(dat)-4
  myp = (numCovs+1)*n.k
  lambda.min.ratio = ifelse(nrow(dat)<myp,0.01,0.0001)
  
  #---------------------------------
  # The sequence of lambda2 : ulam2
  #--------------------------------
  if (is.null(lam2)) {
    lambda_max <- 1- 10^(-3)
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
  
  
  #ulam2[6]<-0.5
  
  initOut = extractMats(dat,n.k=n.k)
  #sparOmega = initOut$sparOmega
  #smoOmega1 = initOut$smoOmega1
  #designMat1 = initOut$designMat1
  #basisMat0 = initOut$basisMat0
  
  
  getSeqLam1HpOut = lapply(as.list(ulam2), function(x){
    getSeqLam1Hp(lambda2=x,dat=dat[,1:2],
                 lambda=lambda, nlam=nlam, sparOmega=initOut$sparOmega, 
                 smoOmega1=initOut$smoOmega1, designMat1=initOut$designMat1,
                 basisMat0=initOut$basisMat0)
  })
  lamGrid= vapply(seq(ulam2), function(i){
    getSeqLam1HpOut[[i]]$ulam
  }, FUN.VALUE =  rep(1, nlam))
  
  
  
  #--------------
  # Step 1, CV fold
  #----------------
  
  dat <- dat[sample(1:nrow(dat), nrow(dat)),]
  foldIndex <-caret::createFolds(dat$Meth_Counts/dat$Total_Counts, k = nfolds)
  
  trainAll <- testAll <- vector("list", nfolds)
  for ( i in seq(nfolds)){
    
    testID <- foldIndex[[i]]
    #trainID <- setdiff(1:nrow(dat), foldIndex[[i]])
    
    trainDat<-dat[-testID,]
    testDat <- dat[testID,]
    
    initOutTrain = extractMats(dat=trainDat,n.k=n.k)
    
    
    # Fit on the train dataset
    
    trainFit <- sparseSmoothGridRaw(dat=trainDat[,1:2], n.k=n.k, ulam2, getSeqLam1HpOut, 
                                    theta=initTheta, stepSize=stepSize, shrinkScale=shrinkScale,
                                    basisMat0=initOutTrain$basisMat0, designMat1=initOutTrain$designMat1, 
                                    numCovs=numCovs, maxInt=maxInt, epsilon=epsilon, accelrt=accelrt, truncation=truncation, mc.cores=mc.cores)
    
    # Calculate the prediction
    
    testPred <- sparseSmoothPred(trainFit=trainFit,trainDatPos=trainDat$Position,testDat=testDat,
                                 basisMat0=initOutTrain$basisMat0,basisMat1=initOutTrain$basisMat1, 
                                 n.k=n.k, numCovs=numCovs,truncation=truncation)
    trainAll[[i]] =trainFit
    testAll[[i]] = testPred
  }
  
  
  # What to export from the cv function
  # best lambda2, best lambda1
  
  #1. calculate the cross validation average/SD matrix of testPred
  
  testPredMean = Reduce("+", testAll) / length(testAll)
  
  testPredSD = apply(simplify2array(testAll), 1:2, sd )
  
  
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
    bestLambda1 = lamGrid[bestInd]
    bestLambda2 = ulam2
  }
  
 
  
  #Hp1 =  (1-bestLambda2)*initOut$sparOmega + bestLambda2*initOut$smoOmega1
  #see1 = getSeqLam1HpOut[[bestInd[2]]]$Hp
  
  initOut = extractMats(dat,n.k=n.k)
  
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
                               lambda2=ulam2[i], Hp=getSeqLam1HpOut[[i]]$Hp, 
                              Linv=getSeqLam1HpOut[[i]]$Linv, 
                               theta=trainAll[[1]]$thetaOut[[i]][,bestIndAllLam2[i]], 
                               stepSize, basisMat0=initOut$basisMat0, designMat1=initOut$designMat1,
                               numCovs, maxInt ,  epsilon , shrinkScale, 
                               accelrt=FALSE, truncation = TRUE, basisMat1 = initOut$basisMat1)
  }, mc.cores=mc.cores)
   
   bestFit <- bestFitAll[[bestInd[2]]]
  
  }
  if(nlam2 == 1){
    bestFit <- sparseSmoothPathRawOneLam1(dat[,1:3], n.k, ulam= bestLambda1, 
                                          lambda2=bestLambda2, Hp=getSeqLam1HpOut[[1]]$Hp, 
                                          Linv=getSeqLam1HpOut[[1]]$Linv, 
                                          theta=trainAll[[1]]$thetaOut[[1]][, bestInd[1]], 
                                          stepSize, basisMat0=initOut$basisMat0, designMat1=initOut$designMat1,
                                          numCovs, maxInt ,  epsilon , shrinkScale, 
                                          accelrt=FALSE, truncation = TRUE, basisMat1 = initOut$basisMat1)
    
  }
  
  if(nlam2==1){
  return(out = list(bestFit = bestFit,testPredMean=testPredMean,testPredSD=testPredSD,
                    ulam2=ulam2, lamGrid=lamGrid, bestLambda1=bestLambda1, bestLambda2=bestLambda2,
                    bestInd=bestInd))
  }else{
    return(out = list(bestFit = bestFit,testPredMean=testPredMean,testPredSD=testPredSD,
                      ulam2=ulam2, lamGrid=lamGrid, bestLambda1=bestLambda1, bestLambda2=bestLambda2,
                      bestInd=bestInd, bestFitAll=bestFitAll)) 
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


