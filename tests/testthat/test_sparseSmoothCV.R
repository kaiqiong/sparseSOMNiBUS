setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat")
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "datSnp5nsig1.RDS", sep = "")
#load("/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")
#saveRDS(dat, path_ref_data  )

dat = readRDS(path_ref_data)

n.snp <- ncol(dat)-4



#library(sparseSOMNiBUS)
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFitPath.R")
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothPred.R")
setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/src")
library(Rcpp)
sourceCpp("sparseOmegaCr.cpp")
sourceCpp("utils.cpp")
sourceCpp("updates.cpp")
sourceCpp("proxGradFit.cpp")

#sourceCpp("sparseSmoothPath.cpp")

n.k = 5
numCovs = ncol(dat)-4

lambda = NULL
nlam = 50

lam2 = NULL
nlam2 = 5

initTheta <- rnorm(n.k*(numCovs+1))
stepSize=0.1
shrinkScale=0.5

maxInt = 1000
epsilon = 1E-6
accelrt = TRUE
truncation = TRUE
mc.cores = 10

nfolds = 10


# Step 1, CV fold



dat <- dat[sample(1:nrow(dat), nrow(dat)),]
foldIndex <-caret::createFolds(dat$Meth_Counts/dat$Total_Counts, k = nfolds)



for ( i in seq(nfolds)){
  
  testID <- foldIndex[[i]]
  #trainID <- setdiff(1:nrow(dat), foldIndex[[i]])
  
  trainDat<-dat[-testID,]
  testDat <- dat[testID,]
  
  initOut = extractMats(dat=trainDat,n.k=n.k)
  

  # Fit on the train dataset
  
  trainFit <- sparseSmoothGrid(dat=trainDat[,1:2], n.k=n.k, lambda=lambda, nlam=nlam, lam2=lam2, nlam2=nlam2, 
                               theta=initTheta, stepSize=stepSize, shrinkScale=shrinkScale,
                   basisMat0=initOut$basisMat0, sparOmega=initOut$sparOmega,
                   smoOmega1=initOut$smoOmega1, designMat1=initOut$designMat1, 
                   numCovs=numCovs, maxInt=maxInt, epsilon=epsilon, accelrt=accelrt, truncation=truncation, mc.cores=mc.cores)
  
  # Calculate the prediction
  
  testPred <- sparseSmoothPred(trainFit=trainFit,trainDatPos=trainDat$Position,testDat=testDat, basisMat0=initOut$basisMat0,
                   basisMat1=initOut$basisMat1, n.k=n.k, numCovs=numCovs,truncation=truncation)

}



i = 1
i = i + 1
points(log(trainFit$lamGrid[,i]),testPred[,i], col = i)

plot(log(trainFit$lamGrid[,1]), testPred[,1], xlim = log(c(min(trainFit$lamGrid), max(trainFit$lamGrid))),
     ylim = c(min(testPred), max(testPred)))
for(i in 2:5){
  points(log(trainFit$lamGrid[,i]), testPred[,i], col = i)
}

trainFit$ulam2
apply(testPred, 2, min)

(lam2Min = which.min(apply(testPred, 2, min)))
# 4
(lamMin = which.min(testPred[,lam2Min])) 
#15


# How about the results using the whole dataset
initOutAll = extractMats(dat,n.k=n.k)



bestFit <- sparseSmoothBest(theta=trainFit$thetaOut[[lam2Min]][, lamMin], stepSize,
                 lambda2=trainFit$ulam2[lam2Min], dat=dat[,1:2], 
                 initOutAll$basisMat0, n.k,
                 initOutAll$sparOmega,initOutAll$smoOmega1,
                 initOutAll$designMat1, lambda =trainFit$lamGrid[lamMin,lam2Min], nlam = 100, numCovs,
                             maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, truncation = TRUE)

fit1 <- fitProxGradCpp(theta, intStepSize = stepSize, lambda1 = trainFit$lamGrid[15,4],
                       dat[,1:2], basisMat0_tilda, n.k,Hp,
                       maxInt, epsilon, shrinkScale,
                       accelrt, numCovs, designMat1_tilda, truncation)

which.min(testPred[,1])
which.min(testPred[,2])
which.min(testPred[,3])
which.min(testPred[,4])
which.min(testPred[,5])

which.min(apply(testPred, 2, min))



min(testPred)
