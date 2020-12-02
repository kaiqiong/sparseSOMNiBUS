setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat")
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "datSnp5nsig1.RDS", sep = "")
#load("/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")
#saveRDS(dat, path_ref_data  )

dat = readRDS(path_ref_data)

n.snp <- ncol(dat)-4



#library(sparseSOMNiBUS)
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFitPath.R")
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

maxInt = 500
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

outlist <- as.list(seq(nfolds))
for( i in seq(length(outlist))){
  outlist[[i]] <-as.list(seq(length(lambda1)*length(lambda2)))
}

gridIndex <- expand.grid(seq(length(lambda1)),seq(length(lambda2)) )

testLossOut <- array(NA, dim= c(4, nfolds, length(lambda1), length(lambda2)),
                     dimnames = list( c("lam1", "lam2", "nzeros", "testLoss"), paste0("Fold", seq(nfolds)),
                                      paste0("lam1", lambda1), paste0("lam2", lambda2)))

# step 1: Spline Basis Set up
# calculate matrices: sparseOmega, smoothOmega1, basisMat0 (intercept), basisMat1 (for rest of covariates)
# These matrices are fixed for fixed n.k and Position
initOut = extractMats(dat=dat,n.k=n.k)

basisMat0 <- initOut$basisMat0
basisMat1 <- initOut$basisMat1
sparOmega <- initOut$sparOmega
smoOmega1 <- initOut$smoOmega1
designMat1 <- initOut$designMat1
#Hp <- (1-lambda2)* sparOmega + lambda2*smoOmega1




stepSize=2
theta_m <- theta <- initTheta <- rep(0, n.k*(ncol(dat)-3))
shrinkScale=0.5

maxInt = 200
epsilon = 1E-6
printDetail = FALSE
accelrt = FALSE

iter = 4

#sourceCpp("utils.cpp")
#sourceCpp("updates.cpp")

theta <- rnorm(length(theta))
truncation= TRUE


# Fit a sequence of lambda1



#sourceCpp("sparseSmoothPath.cpp")

out1 = sparseSmoothGrid(theta, stepSize, lam2 = NULL, nlam2 = 10, dat, basisMat0, n.k, sparOmega,smoOmega1,
                        designMat1, basisMat1,  lambda = NULL, nlam = 100, numCovs,
                        maxInt ,   epsilon , shrinkScale,accelrt=FALSE, 
                        truncation = TRUE, mc.cores = 10)


out1$ulam2


out1$lamGrid[1,]

out1$thetaOut
