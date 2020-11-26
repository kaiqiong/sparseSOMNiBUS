setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat")
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "datSnp5nsig1.RDS", sep = "")
#load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")
#saveRDS(dat, path_ref_data  )

dat = readRDS(path_ref_data)

n.snp <- ncol(dat)-4

# pre-set parameters

n.k = 5
numCovs = ncol(dat)-4
shrinkScale=1/2


lambda2 = 2

lambda1 = 10


source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/utils.R")
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFit.R")
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/oneUpdate.R")
source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothPred.R")
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
Hp <- sparOmega + lambda2*smoOmega1




stepSize=2
theta_m <- theta <- initTheta <- rep(0, n.k*(ncol(dat)-3))
shrinkScale=0.5

maxInt = 200
epsilon = 1E-6
printDetail = FALSE
accelrt = FALSE

iter = 4

#sourceCpp("updates.cpp")

sourceCpp("utils.cpp")

#----------------------------
# test getSeparateThetaCpp()
#---------------------------
theta <- rnorm(length(theta))

see1 = getSeparateThetaCpp(theta, n.k, numCovs)

see2 = getSeparateTheta(theta, n.k, numCovs)



all.equal(see1, see2)


library(microbenchmark)

see = microbenchmark(R = {
  getSeparateTheta(theta, n.k, numCovs)
},
cpp =getSeparateThetaCpp(theta, n.k, numCovs)
, times = 100 )
library(ggplot2)
autoplot(see)
see

#----------------------------
# test estimatePijCpp()
#---------------------------
theta <- rnorm(length(theta))
sourceCpp("utils.cpp")
see1 = estimatePijCpp(theta, basisMat0, designMat1,  n.k, numCovs)

see2 = estimatePij(theta, basisMat0, designMat1, dat,  n.k, numCovs)

all.equal(see1, see2)


all.equal(as.vector(see1$pi_ij), see2$pi_ij)
all.equal(see1$theta.sep, see2$theta.sep)



library(microbenchmark)

see = microbenchmark(R = {
  estimatePij(theta, basisMat0, designMat1, dat,  n.k, numCovs)
},
cpp =estimatePijCpp(theta, basisMat0, designMat1,  n.k, numCovs)
, times = 100 )
library(ggplot2)
autoplot(see)

see
#----------------------------
# test binomObjectCpp()
#---------------------------
theta <- rnorm(length(theta))
sourceCpp("utils.cpp")
truncation = TRUE

basisMatsee = basisMat0

see1 = binomObjectCpp(theta,basisMat0, dat, n.k,numCovs, designMat1, truncation)

all.equal(basisMat0, basisMatsee)

see2 = binomObject(theta,basisMat0, dat, n.k,numCovs, designMat1)

all.equal(see1, see2)


all.equal(as.vector(see1$gNeg2loglik), see2$gNeg2loglik)



library(microbenchmark)

see = microbenchmark(R = {
  binomObject(theta,basisMat0, dat, n.k,numCovs, designMat1)
},
cpp =binomObjectCpp(theta,basisMat0, dat, n.k,numCovs, designMat1, truncation),
times = 100 )

library(ggplot2)
autoplot(see)
see

library(microbenchmark)




see = microbenchmark(R = {
  binomObject(theta,basisMat0, dat, n.k,numCovs, designMat1)
},
cpp =binomObjectCpp(theta,basisMat0, dat, n.k,numCovs, designMat1, truncation),
Naivecpp= binomObjectCppRef(theta,basisMat0, dat, n.k,numCovs, designMat1, truncation),
cppref=binomObjectCpprefineRef(theta,basisMat0, dat, n.k,numCovs, designMat1, truncation),
 times = 100 )

library(ggplot2)
autoplot(see)

#-----------------
# test twoPenaltiesCpp()
#-----------------


Hp <- sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0

iter <- 1

theta <- rnorm(length(theta))
out <- oneUpdate(theta, stepSize, lambda1, dat, basisMat0, n.k, Hp,  numCovs, shrinkScale, 
                 designMat1,theta_m=theta, iter=iter, accelrt)
sourceCpp("utils.cpp")
see1 <- twoPenaltiesCpp(out$theta_l_proximal_sep, lambda1, numCovs, n.k)
see2 <- twoPenalties(out$theta_l_proximal_sep,  lambda1, numCovs, n.k)

all.equal(see1, see2)

library(microbenchmark)
see = microbenchmark(R = {
  twoPenalties(out$theta_l_proximal_sep,  lambda1, numCovs, n.k)
},
cpp =twoPenaltiesCpp(out$theta_l_proximal_sep,  lambda1, numCovs, n.k),
times = 100 )

library(ggplot2)
autoplot(see)
see


#----------------------------
# test lambdaMaxCpp()
#---------------------------

sourceCpp("utils.cpp")
nullFit <- glm(dat$Meth_Counts/dat$Total_Counts~-1+basisMat0,
               family = "binomial", weights = dat$Total_Counts)

mu <- nullFit$fitted.values
see1 = lambdaMax(dat, designMat1, Hp, numCovs, mu)


see2 = lambdaMaxCpp(dat$Meth_Counts,dat$Total_Counts, designMat1, Hp, numCovs, mu)


all.equal(see1, see2)


library(microbenchmark)
see = microbenchmark(R = {
  lambdaMax(dat, designMat1, Hp, numCovs, mu)
},
cpp =lambdaMaxCpp(dat$Meth_Counts,dat$Total_Counts, designMat1, Hp, numCovs, mu),
times = 100 )

library(ggplot2)
autoplot(see)
see

#----------------------------
# test designToTilda()
#---------------------------
sourceCpp("utils.cpp")
L  = chol(Hp)
Hinv = chol2inv(L)
Linv = chol2inv(chol(L))


designMat1 <- c(designMat1, designMat1, designMat1, designMat1)

numCovs <- length(designMat1)
#basisMat0_tilda <- basisMat0 %*% Linv
designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})


see2 = designToTilda(designMat1,Linv,numCovs)
all.equal(designMat1_tilda, see2)

library(microbenchmark)
see = microbenchmark(R = {
  designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})
  
},
cpp =designToTilda(designMat1,Linv,numCovs),
times = 100 )

library(ggplot2)
autoplot(see)
see
# Not much gain by using Cpp -- so stick with R instead


