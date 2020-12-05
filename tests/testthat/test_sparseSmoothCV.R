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
nlam2 = 10

initTheta <- rnorm(n.k*(numCovs+1))
stepSize=0.1
shrinkScale=0.5

maxInt = 10^5
epsilon = 1E-6
accelrt = TRUE
truncation = TRUE
mc.cores = 10

nfolds = 10


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


ulam2

initOut = extractMats(dat,n.k=n.k)
sparOmega = initOut$sparOmega
smoOmega1 = initOut$smoOmega1
designMat1 = initOut$designMat1
basisMat0 = initOut$basisMat0

#'@param ulam2 a sequence of lambda2
getSeqLam1Hp <- function(lambda2,dat, lambda=NULL, nlam, sparOmega, smoOmega1, designMat1, basisMat0){
  Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0
  
  #--- Matrix decomposition for Hp and calculate transformed design matrix
  
  L  = chol(Hp)
  Hinv = chol2inv(L)
  Linv = solve(L)
  #Linv = MASS::ginv(L)

  start_fit <-  getStart(y=dat$Meth_Counts, x=dat$Total_Counts, designMat1 , Hp, Hinv, numCovs, basisMat0)
  
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
  

  return(out = list(Hp=Hp, Linv=Linv, start_fit=start_fit, ulam = ulam, nlam = nlam,
                    L=L, Hinv =Hinv))
    
}

getSeqLam1HpOut = lapply(as.list(ulam2), function(x){
  getSeqLam1Hp(lambda2=x,dat=dat[,1:2],
             lambda=lambda, nlam=nlam, sparOmega=sparOmega, 
             smoOmega1=smoOmega1, designMat1=designMat1,
             basisMat0=basisMat0)
})

# extract the grid values for lambda1
lamGrid= vapply(seq(ulam2), function(i){
  getSeqLam1HpOut[[i]]$ulam
}, FUN.VALUE =  rep(1, nlam))


# Fit a grid of lambdas 

AllOut = parallel::mclapply(seq(ulam2), function(i){
  sparseSmoothPathRaw(dat, n.k, ulam=getSeqLam1HpOut[[i]]$ulam,
                      lambda2=ulam2[i], Hp=getSeqLam1HpOut[[i]]$Hp, 
                      Linv=getSeqLam1HpOut[[i]]$Linv, theta, stepSize, basisMat0, 
                      designMat1,  numCovs, maxInt,  epsilon, shrinkScale,
                      accelrt, truncation)
}, mc.cores=mc.cores)

zeroCovsBool<- thetaOutOri <- thetaOut <-  vector("list", nlam2)
#thetaOut <-  array(NA, c(myp, nlam , nlam2))
lamGrid1 <- matrix(NA, nrow = nlam, ncol = nlam2)
for ( i in seq(ulam2)){
  
  #thetaOut[[i]] <- AllOut[[i]]$thetaMat
  thetaOut[[i]] <- Matrix::Matrix(AllOut[[i]]$thetaMat,sparse = TRUE)
  thetaOutOri [[i]] <- Matrix::Matrix(AllOut[[i]]$thetaMatOriginal,sparse = TRUE)
  # thetaOut[,,i] <- AllOut[[i]]$thetaMat
  lamGrid1[,i] <- AllOut[[i]]$ulam
  zeroCovsBool[[i]] <- AllOut[[i]]$zeroCovsBool
}

length(thetaOut)
all.equal(lamGrid, lamGrid1)

#----- optimcheck for the the Grid fit

checkAlllam <-  vector("list", nlam2)

for ( i in seq(ulam2)){
 
  checkAlllam[[i]]<- 
  vapply(1:nlam, function(j){
    temp = getSeparateThetaCpp(thetaOut[[i]][,j], nk = n.k, numCovs)
    optimcheck(temp, gNeg2loglik=AllOut[[i]]$gNeg2loglik[,j], lambda1=lamGrid[j,i], 
               Hp=getSeqLam1HpOut[[i]]$Hp, 
               L=getSeqLam1HpOut[[i]]$L, Linv=getSeqLam1HpOut[[i]]$Linv, 
               Hpinv=getSeqLam1HpOut[[i]]$Hinv, nk=n.k,eqDelta= 0.01, uneqDelta=10^-5)
  }, FUN.VALUE = rep(TRUE, numCovs))
  
  
}

# End of the check All
#---------------
thetaOut

optimcheck(thetaTildaSep, gNeg2loglik, lambda1, 
           Hp, L, Linv, Hpinv, nk,eqDelta= 10^-5, uneqDelta=10^-5)

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
