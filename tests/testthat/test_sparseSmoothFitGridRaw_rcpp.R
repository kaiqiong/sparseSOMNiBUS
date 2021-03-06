setwd("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat")
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "datSnp5nsig1.RDS", sep = "")
#load("/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")
#saveRDS(dat, path_ref_data  )

dat = readRDS(path_ref_data)

n.snp <- ncol(dat)-4



#library(sparseSOMNiBUS)
source("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFitPath.R")
source("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothPred.R")
source("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/utils.R")
source("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFit.R")
setwd("/scratch/greenwood/kaiqiong.zhao/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/src")
library(Rcpp)
sourceCpp("sparseOmegaCr.cpp")
sourceCpp("utils.cpp")
sourceCpp("updates.cpp")
sourceCpp("proxGradFit.cpp")

#sourceCpp("sparseSmoothPath.cpp")

n.k = 5
numCovs = ncol(dat)-4

lambda = NULL
nlam = 100

lam2 = NULL
nlam2 = 10

initTheta <- rnorm(n.k*(numCovs+1))
stepSize=0.1
shrinkScale=0.5

maxInt = 10^5
epsilon = 1E-6
accelrt = TRUE
truncation = TRUE
mc.cores = 10

nfolds = 5


# Step 0: determine the grid of lambda1 and lambda2 for cross validation for the full dataset

#1. a sequence of lambda2
#2. generate corresponding Hp, L, Hinv, Linv and save them
# save the tilda version of the design matrix as well
#3. calculate lambda.max for lambda1



myp = (numCovs+1)*n.k
lambda.min.ratio = ifelse(nrow(dat)<myp,0.01,0.0001)

#---------------------------------
# The sequence of lambda2 : ulam2
#--------------------------------
if (is.null(lam2)) {
  lambda_max <- 1- 10^(-1)
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


LinvList= lapply(seq(ulam2), function(i){
  getSeqLam1HpOut[[i]]$Linv
})


for ( i in 1:nlam2){
  print((all.equal(LinvList[[i]], getSeqLam1HpOut[[i]]$Linv)))
}

neg2loglikSat = start_fit$neg2loglik_sat
intStepSize = stepSize

sourceCpp("../src_new_big_func/sparseSmoothGridRaw.cpp")
see1 = sparseSmoothGridRawCpp(lamGrid,
                       ulam2,
                       LinvList, theta,  intStepSize,dat, basisMat0,  n.k, maxInt,epsilon,shrinkScale,
                      accelrt, 
                        numCovs, 
                        designMat1, 
                        truncation,
                       neg2loglikSat)





see2 = sparseSmoothGridRaw(dat, n.k, ulam2, getSeqLam1HpOut, theta, stepSize, shrinkScale, basisMat0, designMat1, numCovs,
                    maxInt, epsilon, accelrt, truncation, mc.cores)


all.equal(as.matrix(see2$thetaOut[[1]]), see1$thetaMat[,,1])
all.equal(as.matrix(see2$thetaOut[[2]]), see1$thetaMat[,,2])
all.equal(as.matrix(see2$thetaOut[[3]]), see1$thetaMat[,,3])
all.equal(as.matrix(see2$thetaOut[[4]]), see1$thetaMat[,,4])
all.equal(as.matrix(see2$thetaOut[[5]]), see1$thetaMat[,,5])


all.equal(as.matrix(see2$thetaOutOri[[1]]), as.matrix(see1$thetaMatOriginal[,,1]))
all.equal(as.matrix(see2$thetaOutOri[[2]]), see1$thetaMatOriginal[,,2])
all.equal(as.matrix(see2$thetaOutOri[[3]]), see1$thetaMatOriginal[,,3])
all.equal(as.matrix(see2$thetaOutOri[[4]]), see1$thetaMatOriginal[,,4])
all.equal(as.matrix(see2$thetaOutOri[[5]]), see1$thetaMatOriginal[,,5])










seecomp = microbenchmark(R = {
  sparseSmoothGridRaw(dat, n.k, ulam2, getSeqLam1HpOut, theta, stepSize, shrinkScale, basisMat0, designMat1, numCovs,
                      maxInt, epsilon, accelrt, truncation, mc.cores=1)
    
  }
  ,
cpp ={
  
 
  lamGrid= vapply(seq(ulam2), function(i){
    getSeqLam1HpOut[[i]]$ulam
  }, FUN.VALUE =  rep(1, nlam))
  
  
  LinvList= lapply(seq(ulam2), function(i){
    getSeqLam1HpOut[[i]]$Linv
  })
  
  
  sparseSmoothGridRawCpp(lamGrid,
                            ulam2,
                            LinvList, theta,  intStepSize,dat, basisMat0,  n.k, maxInt,epsilon,shrinkScale,
                            accelrt, 
                            numCovs, 
                            designMat1, 
                            truncation,
                            neg2loglikSat)
}
, times = 5 )
library(ggplot2)
autoplot(seecomp)
seecomp


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


bestInd = which(testPredMean == min(testPredMean), arr.ind = TRUE)

bestLambda1 = lamGrid[bestInd]
bestLambda2 = ulam2[bestInd[2]]





#library(microbenchmark)

#see = microbenchmark(Reduce = {
#  Reduce("+", testAll) / length(testAll)
#},
#apply= apply(simplify2array(testAll), 1:2, mean )

#, times = 10 )
#library(ggplot2)
#autoplot(see)
#see

# Reduce is better than apply

#thetaMatOriginal=thetaMatOriginal,zeroCovsBool=zeroCovsBool

# Extract the estimated penalized coefficients



i = 1
i = i + 1

for( i in 2:100){
  points(uniqPos,  penalBetas[,i], col = i)
}

points(pos, beta.1, col = 'dodgerblue')

#which(trainDat$Position==trainDat$Position[500])


#initOut$basisMat1[500,]==initOut$basisMat1[which(trainDat$Position==trainDat$Position[500])[11],]

length(testAll)

testAll[[1]]


plot(lamGrid)

lamGridPlot = lamGrid
lamGridPlot[1,] = lamGrid[1,]-100000

i = 1
i = i + 1
points(log(lamGridPlot[,i]),testPred[,i], col = i)

plot(log(lamGridPlot[,1]), testPred[,1], xlim = log(c(min(lamGridPlot), max(lamGridPlot))),
     ylim = c(min(testPred), max(testPred)))
for(i in 2:10){
  points(log(lamGridPlot[,i]), testPred[,i], col = i)
}









initOut = extractMats(dat,n.k=n.k)
#Hp1 =  (1-bestLambda2)*initOut$sparOmega + bestLambda2*initOut$smoOmega1
#see1 = getSeqLam1HpOut[[bestInd[2]]]$Hp


bestFit <- sparseSmoothPathRawOneLam1(dat[,1:3], n.k, ulam= bestLambda1, 
                                      lambda2=bestLambda2, Hp=getSeqLam1HpOut[[bestInd[2]]]$Hp, 
                                      Linv=getSeqLam1HpOut[[bestInd[2]]]$Linv, 
                                      theta=trainAll[[1]]$thetaOut[[bestInd[2]]][, bestInd[1]], 
                                      stepSize, basisMat0=initOut$basisMat0, designMat1=initOut$designMat1,
                                      numCovs, maxInt ,  epsilon , shrinkScale, 
                                      accelrt=FALSE, truncation = TRUE, basisMat1 = initOut$basisMat1)

bestFit$zeroCovsBool

plot(bestFit$uniqPos, bestFit$penalBetas)






cvout = sparseSmoothFitCV(dat, n.k, stepSize=0.1, lambda = NULL, nlam = 100, lam2 = NULL, nlam2 = 10, maxInt = 500,
                          epsilon = 1E-6, printDetail = TRUE, initTheta, shrinkScale=0.5,
                          accelrt = TRUE, nfolds = 5, mc.cores)



plotCV(cvout)

cvout$bestFit$zeroCovsBool
