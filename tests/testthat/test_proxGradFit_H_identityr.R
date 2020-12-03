#setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat")
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "datSnp5nsig1.RDS", sep = "")
#load("/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")
#saveRDS(dat, path_ref_data  )

dat = readRDS(path_ref_data)

n.snp <- ncol(dat)-4

# pre-set parameters

n.k = 5
numCovs = ncol(dat)-4
shrinkScale=1/2


lambda2 = 0.2

lambda1 = 10


source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/utils.R")
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFit.R")
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/oneUpdate.R")
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothPred.R")
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFitPath.R")
setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/src")
library(Rcpp)
sourceCpp("sparseOmegaCr.cpp")

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
#theta_m <- theta <- initTheta <- rep(0, n.k*(ncol(dat)-3))
shrinkScale=0.5

maxInt = 200
epsilon = 1E-6
printDetail = FALSE
accelrt = FALSE

iter = 4

sourceCpp("utils.cpp")
sourceCpp("updates.cpp")

theta <- rnorm(length(theta))
truncation= TRUE




#---------------------------------

# test define Hp = identity matrix
#---------------------------------


sourceCpp("proxGradFit.cpp")

truncation= TRUE
theta <- rnorm(length(theta))

lambda2 = 0.5

Hp <- diag(n.k)
lambda1 = 5
seeinit= fitProxGradCpp(theta, stepSize,lambda1=5, dat, basisMat0, n.k, Hp, maxInt = 10^5,
                     epsilon = 1E-20, shrinkScale,accelrt, numCovs, designMat1, basisMat1, truncation)

lambda1 = 3
see2= fitProxGradCpp(seeinit$thetaEst, stepSize,lambda1, dat, basisMat0, n.k, Hp, maxInt = 10^5,
                     epsilon = 1E-20, shrinkScale,accelrt, numCovs, designMat1, basisMat1, truncation)

lambda1 = 1000
see3= fitProxGradCpp(see2$thetaEst, stepSize,lambda1, dat, basisMat0, n.k, Hp, maxInt = 10^5,
                     epsilon = 1E-20, shrinkScale,accelrt, numCovs, designMat1, basisMat1, truncation)
#lambda1 = 4
#see2R= fitProxGrad(seeinit$thetaEst, stepSize,lambda1, dat, basisMat0, n.k, sparOmega, lambda2, smoOmega1, maxInt = 10^3,
#                     epsilon = 1E-20, printDetail = F, shrinkScale,accelrt, numCovs, designMat1, basisMat1)


#Hp <-sparOmega + lambda2*smoOmega1
#see21= fitProxGradCppold(theta, stepSize,lambda1, dat, basisMat0, n.k, sparOmega, lambda2, smoOmega1, maxInt = 10^5,
#                      epsilon = 1E-6, shrinkScale,accelrt, numCovs, designMat1, basisMat1, truncation)



#see2
#all.equal(see2, see21)



ssfit=see2
Hp <- (1-lambda2)* sparOmega + lambda2*smoOmega1
nk = n.k
optimcheck <- function(ssfit, lambda1, Hp, nk,eqDelta= 10^-5, uneqDelta=10^-5){
  
  thetaSep <- ssfit$thetaEstSep[-1]
  zeroCovs <- unlist(lapply( thetaSep, function(x){all(x==0)}))
  indOptimRes <- rep(NA, length(zeroCovs))
  
  zeroCovid <- which(zeroCovs)
  nonzeroCovid <- which(!zeroCovs)
 
  
  for(p in  nonzeroCovid ){
    
   ap = -ssfit$gNeg2loglik[(p*nk+1) :((p+1)*nk)]
    
   # ap =  2* t(designMat1[[p]]) %*% (dat$Meth_Counts-dat$Total_Counts*ssfit$pi_ij) 
    
  top = Hp %*% thetaSep[[p]]
  denom =sqrt(sum(top * thetaSep[[p]] )) 
  term2 = top/denom*lambda1
  
  indOptimRes[p]=sqrt(sum((ap - term2)^2)) <=eqDelta
    
  }
  
  for(p in zeroCovid){
    ap = -ssfit$gNeg2loglik[(p*nk+1) :((p+1)*nk)]
    
    temp1 = Hp %*% ap
   indOptimRes[p]= lambda1 -sqrt(sum(temp1*ap)) >= uneqDelta
  }
  return(indOptimRes)
  
}


optimcheck(see2, lambda1=3, Hp, nk = n.k,eqDelta= 10^-2, uneqDelta=10^-2)


optimcheck(seeinit,lambda1=5, Hp, nk = n.k,eqDelta= 10^-2, uneqDelta=10^-2)

optimcheck(see3, lambda1=1000, Hp, nk = n.k,eqDelta= 10^-2, uneqDelta=10^-10)


intfit = getStart(dat$Meth_Counts, dat$Total_Counts, designMat1, Hp, numCovs,basisMat0)

intfit$lambda_max


#-------------------------------------------------------------
# conclusions: when penalty matrix H is identity matrix, all optimcheck are TRUE !
#-------------------------------------------------------------


