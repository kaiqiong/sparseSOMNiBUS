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

#--- calculate the choleski decomposition and inverse 
L  = chol(Hp)

all.equal(t(L)%*%L, Hp)

Hinv = chol2inv(L)

all.equal(Hinv, solve(Hp))

library(microbenchmark)

see = microbenchmark(solve = {
 solve(Hp)
},
chol ={

  Hinv = chol2inv(chol(Hp))
}
, times = 100 )
library(ggplot2)
autoplot(see)
see

Linv = solve(L)

Linv ==chol2inv(chol(L))

#---- endo of test

L  = chol(Hp)
Hinv = chol2inv(L)
Linv1 = chol2inv(chol(L))
Linv = solve(L)
basisMat0_tilda <- basisMat0 %*% Linv


designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})

all.equal(designMat1_tilda[[1]], designMat1[[1]]%*%Linv)
all.equal(designMat1_tilda[[2]], designMat1[[2]]%*%Linv)
all.equal(designMat1_tilda[[3]], designMat1[[3]]%*%Linv)
all.equal(designMat1_tilda[[4]], designMat1[[4]]%*%Linv)

#seeinit= fitProxGradCpp(theta, stepSize,lambda1, dat, basisMat0, n.k, Hp, maxInt = 10^5,
#                     epsilon = 1E-20, shrinkScale,accelrt, numCovs, designMat1, basisMat1, truncation)

seeinit_correct_design= fitProxGradCpp(theta, stepSize,lambda1, dat, basisMat0_tilda, n.k, Hp, maxInt = 10^5,
                        epsilon = 1E-20, shrinkScale,accelrt, numCovs, designMat1_tilda, Linv, truncation)

optimcheck(seeinit_correct_design,lambda1=lambda1, Hp,L, Linv, nk = n.k, eqDelta = 10^-2,  uneqDelta = 10^-3 )
#optimcheck(seeinit,lambda1=lambda1, Hp,L, Linv, nk = n.k, eqDelta = 10^-2,  uneqDelta = 10^-3 )



lambda1 = 130
see2= fitProxGradCpp(seeinit_correct_design$thetaEstTilda, stepSize,lambda1, dat, basisMat0_tilda, n.k, Hp, maxInt = 10^5,
                     epsilon = 1E-20, shrinkScale,accelrt, numCovs, designMat1_tilda, Linv, truncation)

optimcheck(see2,lambda1=lambda1, Hp,L, Linv, nk = n.k, eqDelta = 10^-2,  uneqDelta = 10^-3 )
#lambda1 = 4
#see2R= fitProxGrad(seeinit$thetaEst, stepSize,lambda1, dat, basisMat0, n.k, sparOmega, lambda2, smoOmega1, maxInt = 10^3,
#                     epsilon = 1E-20, printDetail = F, shrinkScale,accelrt, numCovs, designMat1, basisMat1)


#Hp <-sparOmega + lambda2*smoOmega1
#see21= fitProxGradCppold(theta, stepSize,lambda1, dat, basisMat0, n.k, sparOmega, lambda2, smoOmega1, maxInt = 10^5,
#                      epsilon = 1E-6, shrinkScale,accelrt, numCovs, designMat1, basisMat1, truncation)



#see2
#all.equal(see2, see21)



