setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat")
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
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/lambdaMax.R")
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




#-----------------------------
# test define lambda2 within 0, 1
#------------


sourceCpp("proxGradFit.cpp")

truncation= TRUE
theta <- rnorm(length(theta))



#lambda2 = 1
lambda2 = 0.5
Hp <- (1-lambda2)* sparOmega + lambda2*smoOmega1

startfit = getStart(dat$Meth_Counts, dat$Total_Counts, designMat1, Hp, numCovs, basisMat0)

startfit$lambda_max

lambda1 = 140


# transform the design matrix



L  = chol(Hp)
Hinv = chol2inv(L)
Linv = solve(L)
basisMat0_tilda <- basisMat0 %*% Linv


designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})



seeinit_correct_design= fitProxGradCpp(theta, stepSize,lambda1, dat, basisMat0_tilda, n.k, Hp, maxInt = 10^5,
                        epsilon = 1E-20, shrinkScale,accelrt, numCovs, designMat1_tilda,  truncation)

optimcheck(seeinit_correct_design,lambda1=lambda1, Hp,L, Linv,Hpinv=Hinv, nk = n.k, eqDelta = 10^-2,  uneqDelta = 10^-3 )
#optimcheck(seeinit,lambda1=lambda1, Hp,L, Linv, nk = n.k, eqDelta = 10^-2,  uneqDelta = 10^-3 )



lambda1 = 130
see2= fitProxGradCpp(seeinit_correct_design$thetaEst, stepSize,lambda1, dat, basisMat0_tilda, n.k, Hp, maxInt = 10^5,
                     epsilon = 1E-20, shrinkScale,accelrt, numCovs, designMat1_tilda, truncation)

optimcheck(see2,lambda1=lambda1, Hp,L, Linv, Hinv, nk = n.k, eqDelta = 10^-2,  uneqDelta = 10^-3 )
#lambda1 = 4
#see2R= fitProxGrad(seeinit$thetaEst, stepSize,lambda1, dat, basisMat0, n.k, sparOmega, lambda2, smoOmega1, maxInt = 10^3,
#                     epsilon = 1E-20, printDetail = F, shrinkScale,accelrt, numCovs, designMat1, basisMat1)


#Hp <-sparOmega + lambda2*smoOmega1
#see21= fitProxGradCppold(theta, stepSize,lambda1, dat, basisMat0, n.k, sparOmega, lambda2, smoOmega1, maxInt = 10^5,
#                      epsilon = 1E-6, shrinkScale,accelrt, numCovs, designMat1, basisMat1, truncation)



#see2
#all.equal(see2, see21)



