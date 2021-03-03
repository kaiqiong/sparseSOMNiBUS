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
  lambda_max <- 1- 10^(-5)
  # compute lambda sequence
  #ulam2 <- seq(lambda_max, lambda_max*lambda.min.ratio, length.out = nlam2)
  
  #ulam2 <-  1-exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
  #                   length.out = nlam2))
  
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
  
  initOut = extractMats(dat=trainDat,n.k=n.k)
  
  
  # Fit on the train dataset
  
  trainFit <- sparseSmoothGridRaw(dat=trainDat[,1:2], n.k=n.k, ulam2, getSeqLam1HpOut, 
                               theta=initTheta, stepSize=stepSize, shrinkScale=shrinkScale,
                   basisMat0=initOut$basisMat0, designMat1=initOut$designMat1, 
                   numCovs=numCovs, maxInt=maxInt, epsilon=epsilon, accelrt=accelrt, truncation=truncation, mc.cores=mc.cores)
  
  # Calculate the prediction
  
  testPred <- sparseSmoothPred(trainFit=trainFit,trainDatPos=trainDat$Position,testDat=testDat,
                               basisMat0=initOut$basisMat0,basisMat1=initOut$basisMat1, 
                               n.k=n.k, numCovs=numCovs,truncation=truncation)
  trainAll[[i]] =trainFit
  testAll[[i]] = testPred
}

# Extract the estimated penalized coefficients



dim(trainFit$thetaOutOri[[1]])


uniqPos = unique(trainDat$Position)

uni_rows <- match(uniqPos, trainDat$Position)

dim(initOut$basisMat1[uni_rows,])

penalBetas <- initOut$basisMat1[uni_rows,]%*%trainFit$thetaOutOri[[1]][6:10,]

plot(uniqPos, penalBetas[,1])

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

ulam2
apply(testPred, 2, min)

(lam2Min = which.min(apply(testPred, 2, min)))
# 4
(lamMin = which.min(testPred[,lam2Min])) 
#15


zeroCovsBool[[lam2Min]][,lamMin]

# How about the results using the whole dataset
initOutAll = extractMats(dat,n.k=n.k)





bestFit <- sparseSmoothBest(theta=trainFit$thetaOut[[lam2Min]][, lamMin], stepSize,
                 lambda2=trainFit$ulam2[lam2Min], dat=dat[,1:2], 
                 initOutAll$basisMat0, n.k,
                 initOutAll$sparOmega,initOutAll$smoOmega1,
                 initOutAll$designMat1, lambda =trainFit$lamGrid[lamMin,lam2Min], nlam = 100, numCovs,
                             maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, truncation = TRUE)

bestFit <-fitProxGradCpp(theta=trainFit$thetaOut[[lam2Min]][, lamMin], intStepSize = stepSize,
               lambda1 = lamGrid[lamMin,lam2Min],
               dat[,1:2], basisMat0_tilda, n.k,Hp,
               maxInt, epsilon, shrinkScale,
               accelrt, numCovs, designMat1_tilda, truncation)

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

which.min(testPred[1,])
which.min(testPred[2,])
which.min(testPred[3,])
which.min(testPred[4,])
which.min(testPred[5,])
which.min(testPred[6,])
which.min(testPred[7,])
which.min(testPred[8,])
which.min(testPred[9,])

apply(testPred, 1, which.min)
