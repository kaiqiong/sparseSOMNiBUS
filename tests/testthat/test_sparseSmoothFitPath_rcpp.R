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



#------------------------------------
# Test function sparseSmoothPathCpp()  here the warm start process is implemented in Rcpp instead
#-------------------------------------

sourceCpp("proxGradFit.cpp")
sourceCpp("sparseSmoothPath.cpp")

truncation= TRUE
theta <- rnorm(length(theta))

lambda2 = 0.5



# Call 

see1 = sparseSmoothPath(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                                    designMat1, basisMat1,  lambda = c(180, 198, 200), nlam = 100, numCovs,
                                    maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, 
                        truncation = TRUE, eqDelta=0.01, uneqDelta=10^(-4))



#

sourceCpp("sparseSmoothPath.cpp")
see2= sparseSmoothPathCpp (theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                           designMat1, basisMat1,  lambda = c(0.9,1, 4, 2), nlam = 100, numCovs,
                           maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, truncation = TRUE)
#----- run

lambda = c(180, 198, 200)
myp = (numCovs+1)*n.k
eqDelta=0.01
uneqDelta=10^(-4)


Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0

start_fit <-  getStart(y=dat$Meth_Counts, x=dat$Total_Counts, designMat1 , Hp, numCovs, basisMat0)

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

see = fitProxGradCppSeq(ulam, theta, stepSize, dat, basisMat0, n.k, Hp,
                        maxInt, epsilon, shrinkScale, accelrt, numCovs, designMat1, 
                        basisMat1, truncation,neg2loglikSat=start_fit$neg2loglik_sat, optimcheck)


lambda = c(180, 198, 200)
sourceCpp("sparseSmoothPath.cpp")
see2= sparseSmoothPathCpp (theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                                      designMat1, basisMat1,  lambda = c(0.9,1, 4, 2), nlam = 100, numCovs,
                                      maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, truncation = TRUE)

all.equal(see1$thetaMat, see)




library(microbenchmark)
see = microbenchmark(R = {
  sparseSmoothPath(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                   designMat1, basisMat1,  lambda = c(0.9,1, 4, 2), nlam = 100, numCovs,
                   maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, 
                   truncation = TRUE)
},
cpp ={lambda = c(0.9,1, 4, 2)
myp = (numCovs+1)*n.k
Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0

start_fit <-  getStart(y=dat$Meth_Counts, x=dat$Total_Counts, designMat1 , Hp, numCovs, basisMat0)

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




see = fitProxGradCppSeq(ulam, theta, stepSize, dat, basisMat0, n.k, Hp,
                        maxInt, epsilon, shrinkScale, accelrt, numCovs, designMat1, 
                        basisMat1, truncation,neg2loglikSat=start_fit$neg2loglik_sat, optimcheck)},
times = 10 )

library(ggplot2)
autoplot(see)
see
