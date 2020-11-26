


#@param y meth_counts
#@param x total_counts

getStart <- function(y, x, designMat1, Hp, numCovs, basisMat0){
  
  nullFit <- glm(cbind(y, x-y)~-1+basisMat0, family = "binomial")
  
  mu <- nullFit$fitted.values
  
  lambda_max <-   lambdaMax(y, x, designMat1, Hp, numCovs, mu)
  
  nulldev <- nullFit$deviance
  
  #nega2lik <- (-2)*sum(dat$Meth_Counts*log(mu)+(dat$Total_Counts-dat$Meth_Counts)*log(1-mu))
  eps <- 10 * .Machine$double.eps
  mu_s <- y/x
  if(truncation){
    mu_s[which(mu_s > 1 - eps)] <- 1 - eps
    mu_s[which(mu_s < eps)] <- eps
  }
  neg2loglik_sat <- (-2)*sum(y*log(mu_s)+(x-dat$Meth_Counts)*log(1-mu_s))
  
  list(nulldev = nulldev, mu= mu, lambda_max = lambda_max, neg2loglik_sat=neg2loglik_sat)
  
}


#lambda2 = 0.5
#Hp = (1-lambda2)*sparOmega + lambda2*smoOmega1
#lambda2 = 200
#Hp = sparOmega + lambda2*smoOmega1
#getStart(dat$Meth_Counts, dat$Total_Counts, designMat1, Hp, numCovs, basisMat0 )$lambda_max

#fitProxGradCpp

#if (is.null(lambda)) {
#  if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")

# compute lambda max: to add code here
#  lambda_max <- start_val$lambda_max

# compute lambda sequence
#  ulam <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
#                  length.out = nlam))
#} else { # user provided lambda values
#  user_lambda = TRUE
#  if (any(lambda < 0)) stop("lambdas should be non-negative")
#  ulam = as.double(rev(sort(lambda)))
#  nlam = as.integer(length(lambda))
#}



sparseSmoothPath <- function(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                             designMat1, basisMat1,  lambda = NULL, nlam = 100, numCovs,
                             maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, truncation = TRUE,
                             eqDelta, uneqDelta){
  myp = (numCovs+1)*n.k
  
  Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0
  
  start_fit <-  getStart(y=dat$Meth_Counts, x=dat$Total_Counts, designMat1 , Hp, numCovs, basisMat0)
  
  
  #--- Matrix decomposition for Hp and calculate transformed design matrix
  
  L  = chol(Hp)
  Hinv = chol2inv(L)
  Linv = solve(L)
  
  basisMat0_tilda <- basisMat0 %*% Linv # this tilde does depend on L-- H -- lambda2, should be calculated for each lambda2
  designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})
  
  #-------------------
  
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
  
  lossVec <- rep(NA,nlam)
  thetaMat <- matrix(NA, nrow = nlam, ncol = myp )
  checkall <- matrix(NA, nrow = nlam, ncol = numCovs)
  fit1 <- fitProxGradCpp(theta, intStepSize = stepSize, lambda1 = ulam[1]+100000, dat, basisMat0_tilda, n.k,Hp,
                         maxInt, epsilon, shrinkScale,
                         accelrt, numCovs, designMat1_tilda, truncation)
  
  
  lossVec[1] <- fit1$lossSum
  thetaMat[1,] <- fit1$thetaEst
  
  checkall[1,] <-  optimcheck(fit1, ulam[1], Hp, L, Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
  for( i in 2:length(ulam)){
    fit1 <- fitProxGradCpp(fit1$thetaEst, intStepSize = stepSize, lambda1 = ulam[i], dat, basisMat0_tilda, n.k,Hp,
                           maxInt, epsilon, shrinkScale,
                           accelrt, numCovs, designMat1_tilda, truncation)
    lossVec[i] <- fit1$lossSum
    thetaMat[i,] <- fit1$thetaEst
    
    checkall[i,] <-  optimcheck(fit1, ulam[1], Hp, L,Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
    
    
    if(fit1$neg2loglik < start_fit$neg2loglik_sat) break
    
  }
  
  return(out = list(thetaMat=thetaMat, lossVec=lossVec, ulam= ulam, checkall = checkall))
  
  
}



sparseSmoothPathCpp <- function(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                                designMat1, basisMat1,  lambda = NULL, nlam = 100, numCovs,
                                maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, truncation = TRUE){
  myp = (numCovs+1)*n.k
  Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0
  
  start_fit <-  getStart(y=dat$Meth_Counts, x=dat$Total_Counts, designMat1 , Hp, numCovs, basisMat0)
  
  #--- Matrix decomposition for Hp and calculate transformed design matrix
  
  L  = chol(Hp)
  Hinv = chol2inv(L)
  Linv = solve(L)
  
  basisMat0_tilda <- basisMat0 %*% Linv # this tilde does depend on L-- H -- lambda2, should be calculated for each lambda2
  designMat1_tilda <- lapply(designMat1, function(x){x%*%Linv})
  
  
  
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
  
  see = fitProxGradCppSeq(ulam, theta, stepSize, dat, basisMat0_tilda, n.k, Hp,
                          maxInt, epsilon, shrinkScale, accelrt, numCovs, designMat1_tilda, 
                          truncation,neg2loglikSat=start_fit$neg2loglik_sat, Linv)
  
  
  return(see)
  
  
}


