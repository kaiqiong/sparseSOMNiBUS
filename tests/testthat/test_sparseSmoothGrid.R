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


library(sparseSOMNiBUS)
setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/src")
library(Rcpp)
sourceCpp("sparseOmegaCr.cpp")

sourceCpp("utils.cpp")

sourceCpp("updates.cpp")

sourceCpp("proxGradFit.cpp")
sourceCpp("sparseSmoothPath.cpp")

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
