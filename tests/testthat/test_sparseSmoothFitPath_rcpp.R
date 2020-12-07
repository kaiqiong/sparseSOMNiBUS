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



#------------------------------------
# Test function sparseSmoothPathCpp()  here the warm start process is implemented in Rcpp instead
#-------------------------------------

#sourceCpp("proxGradFit.cpp")
#sourceCpp("sparseSmoothPath.cpp")

truncation= TRUE
theta <- rnorm(length(theta))

lambda2 = 0.5



# Call 

time0 = Sys.time()

lambda = c(180, 198, 200)
myp = (numCovs+1)*n.k
eqDelta=0.01
uneqDelta=10^(-4)

Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0



#--- Matrix decomposition for Hp and calculate transformed design matrix

L  = chol(Hp)
Hinv = chol2inv(L)
Linv = solve(L)



basisMat0_tilda <- basisMat0 %*% Linv # this tilde does depend on L-- H -- lambda2, should be calculated for each lambda2
designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})

#-------------------
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

# Fit a sequence of 

#lossVec <- rep(NA,nlam)

thetaMat <- matrix(NA, nrow = myp, ncol = nlam )
# checkall <- matrix(NA, nrow = numCovs, ncol =nlam)



fit1 <- fitProxGradCpp(theta, intStepSize = stepSize, lambda1 = ulam[1]+100000, dat, basisMat0_tilda, n.k,Hp,
                       maxInt, epsilon, shrinkScale,
                       accelrt, numCovs, designMat1_tilda, truncation)


#lossVec[1] <- fit1$lossSum
thetaMat[,1] <- fit1$thetaEst

# checkall[,1] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[1], Hp, L, Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
for( i in 2:length(ulam)){
  fit1 <- fitProxGradCpp(fit1$thetaEst, intStepSize = stepSize, lambda1 = ulam[i], dat, basisMat0_tilda, n.k,Hp,
                         maxInt, epsilon, shrinkScale,
                         accelrt, numCovs, designMat1_tilda, truncation)
  #lossVec[i] <- fit1$lossSum
  thetaMat[,i] <- fit1$thetaEst
  
  #  checkall[,i] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[i], Hp, L,Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
  
  
  if(fit1$neg2loglik < start_fit$neg2loglik_sat) break
  
}



ulam = c(ulam, 170, 160, 150)
nlam = length(ulam)
library(microbenchmark)

see = microbenchmark(R = {
  thetaMat <- matrix(NA, nrow = myp, ncol = nlam )
 # checkall <- matrix(NA, nrow = numCovs, ncol =nlam)
  fit1 <- fitProxGradCpp(theta, intStepSize = stepSize, lambda1 = ulam[1]+100000, dat, basisMat0_tilda, n.k,Hp,
                         maxInt, epsilon, shrinkScale,
                         accelrt, numCovs, designMat1_tilda, truncation)
  
  
  #lossVec[1] <- fit1$lossSum
  thetaMat[,1] <- fit1$thetaEst
  
 # checkall[,1] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[1], Hp, L, Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
  for( i in 2:length(ulam)){
    fit1 <- fitProxGradCpp(fit1$thetaEst, intStepSize = stepSize, lambda1 = ulam[i], dat, basisMat0_tilda, n.k,Hp,
                           maxInt, epsilon, shrinkScale,
                           accelrt, numCovs, designMat1_tilda, truncation)
    #lossVec[i] <- fit1$lossSum
    thetaMat[,i] <- fit1$thetaEst
    
  #  checkall[,i] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[i], Hp, L,Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
    
    
    if(fit1$neg2loglik < start_fit$neg2loglik_sat) break
    
  }
  
},
cpp =fitProxGradCppSeq(ulam, theta, stepSize, dat, basisMat0_tilda, n.k, Hp,
                       maxInt, epsilon, shrinkScale, accelrt, numCovs, designMat1_tilda, 
                       truncation,neg2loglikSat=start_fit$neg2loglik_sat)
, times = 10 )
library(ggplot2)
autoplot(see)
see


# conclusiong:  write an fitProxGradCpp Seq doesn't seem to add a lot speed benefit
# similar computational time for 3 lambda
#-----------------
time0= Sys.time()
see1 = sparseSmoothPath(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                                    designMat1, basisMat1,  lambda = c( 200,300), nlam = 100, numCovs,
                                    maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, 
                        truncation = TRUE, eqDelta=0.01, uneqDelta=10^(-4))
print(Sys.time()-time0)


time0= Sys.time()
see1 = sparseSmoothPath(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                        designMat1, basisMat1,  lambda = NULL, nlam = 100, numCovs,
                        maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, 
                        truncation = TRUE)
print(Sys.time()-time0)
#Time difference of 7.464239 mins

time0= Sys.time()
see1 = sparseSmoothPath(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                        designMat1, basisMat1,  lambda = NULL, nlam = 100, numCovs,
                        maxInt = 200,  epsilon = 1E-6, shrinkScale,accelrt=FALSE, 
                        truncation = TRUE)
print(Sys.time()-time0)
#Time difference of 15.15335 secs
time0= Sys.time()
see1 = sparseSmoothPath(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                        designMat1, basisMat1,  lambda = NULL, nlam = 100, numCovs,
                        maxInt = 10^5,  epsilon = 1E-6, shrinkScale,accelrt=FALSE, 
                        truncation = TRUE)
print(Sys.time()-time0)
#Time difference of 15.60009 secs
#---Another case ---
lambda = NULL
nlam = 100

Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0



#--- Matrix decomposition for Hp and calculate transformed design matrix

L  = chol(Hp)
Hinv = chol2inv(L)
Linv = solve(L)



basisMat0_tilda <- basisMat0 %*% Linv # this tilde does depend on L-- H -- lambda2, should be calculated for each lambda2
designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})

#-------------------
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

# Fit a sequence of fitProxGradCppSeq
time0 = Sys.time()
res_now = fitProxGradCppSeq(ulam, theta, stepSize, dat, basisMat0_tilda, n.k, Hp,
                        maxInt, epsilon, shrinkScale, accelrt, numCovs, designMat1_tilda, 
                        truncation,neg2loglikSat=start_fit$neg2loglik_sat)


res_now1 = fitProxGradCppSeq(ulam, theta, stepSize, dat, basisMat0_tilda, n.k, Hp,
                            maxInt, epsilon, shrinkScale, accelrt, numCovs, designMat1_tilda, 
                            truncation,neg2loglikSat=start_fit$neg2loglik_sat)


all.equal(res_now1$thetaMat, t(res_now$thetaMat) )


all.equal(res_now1$gNeg2loglikTilda, t(res_now$gNeg2loglikTilda) )

print(Sys.time()-time0)
#Time difference of 5.576715 mins

dim(res_now$thetaMat)

checkAlllam <- matrix(NA, nrow = nlam, ncol = numCovs)

lamv <- ulam
lamv[1]<- ulam[1]+100000
for ( i in seq(length(ulam))){
  thetaSep =  getSeparateThetaCpp(res_now$thetaMat[,i], n.k, numCovs)
  
  checkAlllam[i,] <- optimcheck(thetaSep, res_now$gNeg2loglikTilda[,i], lamv[i], Hp, L, Linv, Hinv, n.k, eqDelta = 0.1, uneqDelta = 10^-5 )
}


checkAlllam





#-------------------------
# compare directly run codes v.s wrap code within a R function ""
#-------------------------


lambda = c(400, 300, 200, 100)
see = microbenchmark(R = {
  Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0
  #--- Matrix decomposition for Hp and calculate transformed design matrix
  
  L  = chol(Hp)
  Hinv = chol2inv(L)
  Linv = solve(L)

  
  basisMat0_tilda <- basisMat0 %*% Linv # this tilde does depend on L-- H -- lambda2, should be calculated for each lambda2
  designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})
  
  #-------------------
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
  
  # Fit a sequence of 
  
  #lossVec <- rep(NA,nlam)
  thetaMat <- matrix(NA, nrow = myp, ncol = nlam )
  checkall <- matrix(NA, nrow = numCovs, ncol =nlam)
  
  fit1 <- fitProxGradCpp(theta, intStepSize = stepSize, lambda1 = ulam[1]+100000, dat, basisMat0_tilda, n.k,Hp,
                         maxInt, epsilon, shrinkScale,
                         accelrt, numCovs, designMat1_tilda, truncation)
  
  
  #lossVec[1] <- fit1$lossSum
  thetaMat[,1] <- fit1$thetaEst
  
  checkall[,1] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[1], Hp, L, Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
  for( i in 2:length(ulam)){
    fit1 <- fitProxGradCpp(fit1$thetaEst, intStepSize = stepSize, lambda1 = ulam[i], dat, basisMat0_tilda, n.k,Hp,
                           maxInt, epsilon, shrinkScale,
                           accelrt, numCovs, designMat1_tilda, truncation)
    #lossVec[i] <- fit1$lossSum
    thetaMat[,i] <- fit1$thetaEst
    
    checkall[,i] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[i], Hp, L,Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
    
    
    if(fit1$neg2loglik < start_fit$neg2loglik_sat) break
    
  }
  
},
Rfun =sparseSmoothPath(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                       designMat1, basisMat1,  lambda , nlam = 100, numCovs,
                       maxInt,  epsilon, shrinkScale,accelrt=FALSE, 
                       truncation = TRUE, eqDelta=0.01, uneqDelta=10^(-4))
, times = 5 )
library(ggplot2)
autoplot(see)
see


# The comparison is using different values of epsilon and maxInt make them comparable

# True conclusion: No differences
 
#--------
# Compare .Call and directly fitProxGrad

# ---- no difference in terms of computational time
lambda = c(400, 300, 200, 100)
see = microbenchmark(functionName = {
  sparseSmoothPathOld(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                   designMat1, basisMat1,  lambda , nlam = 100, numCovs,
                   maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, 
                   truncation = TRUE, eqDelta=0.01, uneqDelta=10^(-4))
  
},
Call =sparseSmoothPath(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                       designMat1, basisMat1,  lambda , nlam = 100, numCovs,
                       maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, 
                       truncation = TRUE, eqDelta=0.01, uneqDelta=10^(-4))
, times = 5 )
library(ggplot2)
autoplot(see)
see


print(Sys.time()-time0)

#

sourceCpp("sparseSmoothPath.cpp")
see2= sparseSmoothPathCpp (theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                           designMat1, basisMat1,  lambda = c(0.9,1, 4, 2), nlam = 100, numCovs,
                           maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, truncation = TRUE)

