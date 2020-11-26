

dat <- readRDS("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat/data/datSnp5nsig1.RDS")

setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/src")
library(Rcpp)
sourceCpp("sparseOmegaCr.cpp")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/utils.R")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFit.R")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/oneUpdate.R")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothPred.R")


set.seed(1231)
# randomly shuffle the rows of dat matrix
dat <- dat[sample(1:nrow(dat), 100),]


numCovs = ncol(dat)-4
shrinkScale=1/2

n.k = 5
lambda1 = 0.1

lambda2 = 0.1
stepSize=2
theta <- initTheta <- rep(0, n.k*(ncol(dat)-3))
shrinkScale=0.5

maxInt = 200
epsilon = 1E-6
printDetail = FALSE
accelrt = FALSE

lambda1= seq(0, 2, 0.1)
lambda2=seq(0, 2, 0.1)

#----------
library(foreach)
library(doParallel)
numCores=detectCores()
registerDoParallel(2)
#-------

setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Simu")
sparseSmoothFitCV <- function(dat, n.k, stepSize=0.1,lambda1= seq(0, 5, 0.1), lambda2=seq(0, 5, 0.1), maxInt = 500,
                              epsilon = 1E-6, printDetail = TRUE, initTheta, shrinkScale=0.5,
                              accelrt = TRUE, nfolds = 10){
}
  
  foldIndex <-caret::createFolds(dat$Meth_Counts/dat$Total_Counts, k = 10)
  
  total_index <- 1:length(foldIndex)
  sep_index <- split(total_index, ceiling(total_index/2))
  
  for(jj in 1:length(sep_index)){
    RES <- foreach( i = sep_index[[jj]], .combine=rbind) %dopar% {
      testID <- foldIndex[[i]]
      trainID <- setdiff(1:nrow(dat), foldIndex[[i]])
      
      trainDat<-dat[trainID,]
      testDat <- dat[testID,]
      for ( lam1 in lambda1){
        for(lam2 in lambda2){
          
          fitnow <- sparseSmoothFit(dat=trainDat, n.k=n.k, stepSize=stepSize,lambda1=lam1, lambda2=lam2,maxInt=maxInt,epsilon=epsilon,
                                    printDetail=printDetail, initTheta=initTheta,shrinkScale=shrinkScale)
          # test on the test dataset
          
          testMats <- extractMats(dat=testDat, n.k=n.k)
          numCovs <- testMats$numCovs
          Hp <- testMats$sparOmega + lam2*testMats$smoOmega1
          testOut <- binomObject(theta=fitnow$thetaEst,basisMat0=testMats$basisMat0,dat=testDat,n.k=n.k,
                                 numCovs=numCovs,designMat1=testMats$designMat1)
          penTerms <- twoPenalties( getSeparateTheta(fitnow$thetaEst, n.k = n.k, numCovs = numCovs),
                                    Hp=Hp, lambda1=lam1, numCovs=testMats$numCovs, n.k=n.k)
          
     out <-  c( fold = i, lam1=lam1, lam2=lam2, testLoss=testOut$neg2loglik + penTerms, binomLoss = testOut$neg2loglik,
                           penTerms=penTerms)
     save(out, fitnow,  file = paste0("Res_fold_", i, "_lam1_", match(lam1, lambda1), "_lam2_", match(lam2, lambda2), ".RDATA"))
          
          
    }
  }


      }
    }
    
  
 # see = readRDS (file = paste0("Res_fold_", i, "_lam1_", match(lam1, lambda1), "_lam2_", match(lam2, lambda2), ".RDA") )
  