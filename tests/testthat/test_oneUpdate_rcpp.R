setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat")
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "datSnp5nsig1.RDS", sep = "")
#load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")
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
lambda2 = 0.5
Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1




stepSize=2
theta_m <- theta <- initTheta <- rep(0, n.k*(ncol(dat)-3))
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



#------------------------------------
# Test function proximalOperatorCpp()
#-------------------------------------


t=0.0001
lambda1=0.0002

t = 2000
lambda1 = 1000
u_p=theta[1:n.k]
theta <- rnorm(length(theta))
u_p=theta[1:n.k]
sum( (u_p %*% Hp) * u_p)

x = u_p%*% Hp
y = u_p

library(microbenchmark)
see = microbenchmark(R = {
  sum(x *y)
},
cpp_std = dot_std(x,y),
cpp_arma = dot_arma(x,y)
, times = 100 )
library(ggplot2)
autoplot(see)

see

sourceCpp("updates.cpp")




see1 = proximalOperatorCpp(t, lambda1, u_p, n.k)
see2 = proximalOperator(t, lambda1, u_p)

see1a = proximalOperatorCppOld(t, lambda1, Hp, u_p, n.k)
all.equal(see1, see2)

all.equal(as.vector(see1), see2)


library(microbenchmark)

see = microbenchmark(R = {
  proximalOperator(t, lambda1, Hp, u_p)
},
cpp =proximalOperatorCpp(t, lambda1, Hp, u_p, n.k),
cpp_ARMA = proximalOperatorCppOld(t, lambda1, Hp, u_p, n.k)
, times = 100 )
library(ggplot2)
autoplot(see)
see

#------------------------------------
# Test function thetaUpdateCpp()
#-------------------------------------
sourceCpp("updates.cpp")

binom_out <- binomObject(theta,basisMat0, dat, n.k,numCovs, designMat1)
gBinomLoss <- binom_out$gNeg2loglik
iter
truncation = T
see1 = thetaUpdateCpp(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, dat, designMat1, 
                      theta_m=theta, iter,accelrt, truncation)
see2 = thetaUpdate(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, dat, designMat1, 
                   theta_m=theta, iter,accelrt, truncation)


all.equal(see1,see2)

all.equal(as.vector(see1$binom_out_new$gNeg2loglik), see2$binom_out_new$gNeg2loglik)
#expect_true(all.equal(see1,see2))

library(microbenchmark)

see = microbenchmark(R = {
  thetaUpdate(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, dat, designMat1, 
              theta_m=theta, iter,accelrt, truncation)
 
},
cpp = thetaUpdateCpp(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, dat, designMat1, 
                           theta_m=theta, iter,accelrt, truncation)

, times = 100 )
library(ggplot2)
autoplot(see)
see

#----------------------------------------
# Test function oneUpdateCpp
#----------------------------------------

theta <- rnorm(length(theta))
truncation= TRUE

see1 = oneUpdateCpp(theta, stepSize=10, lambda1, dat, basisMat0, nk=n.k, Hp, numCovs, shrinkScale, designMat1,
                    theta, iter = 4 , accelrt, truncation)

see2=oneUpdate(theta, stepSize=10, lambda1, dat, basisMat0, n.k, Hp, numCovs,
               shrinkScale,designMat1, theta, iter, accelrt)

test_that("Rcpp returns the save backtracking line search output as R function ", {
  expect_equal(all.equal(as.vector(see1$G_t_theta), see2$G_t_theta))
  expect_equal(all.equal(as.vector(see1$theta_l_proximal), see2$theta_l_proximal))
})
all.equal(see1, see2)

all.equal(as.vector(see1$G_t_theta), see2$G_t_theta)
all.equal(as.vector(see1$theta_l_proximal), see2$theta_l_proximal)

library(microbenchmark)

see = microbenchmark(R = {
  oneUpdate(theta, stepSize=2, lambda1, dat, basisMat0, n.k, Hp, numCovs,
            shrinkScale,designMat1, theta_m, iter, accelrt)
},
Rcpp =oneUpdateCpp(theta, stepSize=2, lambda1, dat, basisMat0, nk=n.k, Hp, numCovs, shrinkScale, designMat1,
                  theta_m, iter = 4 , accelrt,  truncation)
, times = 100 )
library(ggplot2)

pdf("Computational_time_for_one_ProxGrad_update.pdf", width=8, height = 5)
autoplot(see)+ggtitle("Computation time for one update/step of the proximal gradient descent algorithm")
dev.off()
see

  