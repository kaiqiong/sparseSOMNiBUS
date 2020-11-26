setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat")
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


source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/utils.R")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFit.R")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/oneUpdate.R")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothPred.R")

setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/src")
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

sourceCpp("updates.cpp")

see1 = oneUpdateCpp(theta, stepSize=10, lambda1, dat, basisMat0, nk=n.k, Hp, numCovs, shrinkScale, designMat1,
             theta_m, iter = 4 , accelrt, binomObject, thetaUpdate)

see2=oneUpdate(theta, stepSize=10, lambda1, dat, basisMat0, n.k, Hp, numCovs,
                shrinkScale,designMat1, theta_m, iter, accelrt)

test_that("Rcpp returns the save backtracking line search output as R function ", {
  expect_equal(all.equal(as.vector(see1$G_t_theta), see2$G_t_theta))
  expect_equal(all.equal(as.vector(see1$theta_l_proximal), see2$theta_l_proximal))
})
all.equal(see1, see2)

all.equal(as.vector(see1$G_t_theta), see2$G_t_theta)
all.equal(as.vector(see1$theta_l_proximal), see2$theta_l_proximal)

library(microbenchmark)

see = microbenchmark(R = {
  oneUpdate(theta, stepSize=200, lambda1, dat, basisMat0, n.k, Hp, numCovs,
            shrinkScale,designMat1, theta_m, iter, accelrt)
},
                     cpp =oneUpdateCpp(theta, stepSize=200, lambda1, dat, basisMat0, nk=n.k, Hp, numCovs, shrinkScale, designMat1,
                                       theta_m, iter = 4 , accelrt, binomObject, thetaUpdate)
                     , times = 100 )
library(ggplot2)
autoplot(see)

