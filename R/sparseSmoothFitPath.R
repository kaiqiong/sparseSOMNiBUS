#' @useDynLib sparseSOMNiBUS
#' @importFrom Rcpp sourceCpp
#> NULL


#@param y meth_counts
#@param x total_counts

getStart <- function(y, x, designMat1, Hp,Hp_inv, numCovs, basisMat0){
  
  nullFit <- glm(cbind(y, x-y)~-1+basisMat0, family = "binomial")
  
  mu <- nullFit$fitted.values
  
  lambda_max <-   lambdaMax(y, x, designMat1, Hp,Hp_inv, numCovs, mu)
  
  nulldev <- nullFit$deviance
  
  #nega2lik <- (-2)*sum(dat$Meth_Counts*log(mu)+(dat$Total_Counts-dat$Meth_Counts)*log(1-mu))
  eps <- 10 * .Machine$double.eps
  mu_s <- y/x
  if(truncation){
    mu_s[which(mu_s > 1 - eps)] <- 1 - eps
    mu_s[which(mu_s < eps)] <- eps
  }
  neg2loglik_sat <- (-2)*sum(y*log(mu_s)+(x-y)*log(1-mu_s))
  
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
                             designMat1, lambda = NULL, nlam = 100, numCovs,
                             maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, truncation = TRUE){
  
  Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0
  
  #--- Matrix decomposition for Hp and calculate transformed design matrix
  
  L  = chol(Hp)
  Hinv = chol2inv(L)
  Linv = solve(L)
  #Linv = MASS::ginv(L)
  
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
  gNeg2loglikTilda <- matrix(NA, nrow = myp, ncol =nlam);
  lamnow = ulam[1]+100000
  
  fit1 <- fitProxGradCpp(theta, intStepSize = stepSize, lambda1 = ulam[1]+100000, dat[,1:2], basisMat0_tilda, n.k,Hp,
                         maxInt, epsilon, shrinkScale,
                         accelrt, numCovs, designMat1_tilda, truncation)
  
  
  #lossVec[1] <- fit1$lossSum
  thetaMat[,1] <- fit1$thetaEst
  gNeg2loglikTilda[,1]  <-  fit1$gNeg2loglik
  
  #checkall[,1] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[1], Hp, L, Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
  for( i in 2:length(ulam)){
    
    fit1 <- fitProxGradCpp(fit1$thetaEst, intStepSize = stepSize, lambda1 = ulam[i], dat[,1:2], basisMat0_tilda, n.k,Hp,
                           maxInt, epsilon, shrinkScale,
                         accelrt, numCovs, designMat1_tilda, truncation)
    thetaMat[,i] <- fit1$thetaEst
    gNeg2loglikTilda[,i]  <-  fit1$gNeg2loglik
   # checkall[,i] <-  optimcheck(fit1$thetaEstSep, fit1$gNeg2loglik, ulam[i], Hp, L,Linv, Hpinv = Hinv, n.k, eqDelta, uneqDelta )
    
    
    if(fit1$neg2loglik < start_fit$neg2loglik_sat) break
    
  }
  thetaMatOriginal <- 
  vapply(seq(ulam), function(i){
    unlist(lapply(getSeparateThetaCpp(thetaMat[,i], n.k, numCovs), function(x){Linv%*%x}))
  }, FUN.VALUE = rep(1, myp))
  
  
  
  return(out = list(thetaMat=thetaMat, ulam= ulam, gNeg2loglik=gNeg2loglikTilda, thetaMatOriginal=thetaMatOriginal))
  
}

#'@param dat:  columns of outcomes: "Meth_Counts" and "Total_Counts"
sparseSmoothGrid <- function(dat, n.k, lambda = NULL, nlam = 100, lam2 = NULL,  nlam2= 10, theta, stepSize, shrinkScale,
                             basisMat0, sparOmega,smoOmega1, designMat1, numCovs,
                             maxInt = 10^5,  epsilon = 1E-20, accelrt=FALSE, truncation = TRUE,
                             mc.cores = 10){
  myp = (numCovs+1)*n.k
  lambda.min.ratio = ifelse(nrow(dat)<myp,0.01,0.0001)
  
  #---------------------------------
  # The sequence of lambda2 : ulam2
  #--------------------------------
  if (is.null(lam2)) {
    lambda_max <- 0.999
    # compute lambda sequence
    ulam2 <- seq(lambda_max, lambda_max*lambda.min.ratio, length.out = nlam2)
      #exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
      #             length.out = nlam2-1))
    #ulam2 <- c(ulam2, 0)
  } else { # user provided lambda values
    user_lambda2 = TRUE
    if (any(lam2 < 0)) stop("lambdas should be non-negative")
    ulam2 = as.double(rev(sort(lam2)))
    nlam2 = as.integer(length(lam2))
  }
  #---------------------------------
  # The length of the sequence of lambda1: nlam or the length of provided lambda
  #--------------------------------
  if(!is.null(lambda)){
    nlam = as.integer(length(lambda))
  }
  
  thetaOutOri <- thetaOut <-  vector("list", nlam2)
  #thetaOut <-  array(NA, c(myp, nlam , nlam2))
  lamGrid <- matrix(NA, nrow = nlam, ncol = nlam2)
  #dimnames(lamGrid) <- list("lambda", "alpha")
  
  AllOut = parallel::mclapply(seq(ulam2), function(i){
    sparseSmoothPath(theta, stepSize,lambda2=ulam2[i], dat, basisMat0, n.k, sparOmega,smoOmega1,
                            designMat1,   lambda, nlam, numCovs,
                            maxInt ,  epsilon , shrinkScale,accelrt, 
                            truncation)

  }, mc.cores=mc.cores)
  
  for ( i in seq(ulam2)){
    
    #thetaOut[[i]] <- AllOut[[i]]$thetaMat
    thetaOut[[i]] <- Matrix::Matrix(AllOut[[i]]$thetaMat,sparse = TRUE)
    thetaOutOri [[i]] <- Matrix::Matrix(AllOut[[i]]$thetaMatOriginal,sparse = TRUE)
    # thetaOut[,,i] <- AllOut[[i]]$thetaMat
    lamGrid[,i] <- AllOut[[i]]$ulam
  }
  
  # checkall <- matrix(NA, nrow = numCovs, ncol =nlam)
  
  
  return(out = list(thetaOut=thetaOut, lamGrid= lamGrid, ulam2 = ulam2, thetaOutOri=thetaOutOri))
 
}





sparseSmoothPathCpp <- function(theta, stepSize, lambda2=0.5, dat, basisMat0, n.k, sparOmega,smoOmega1,
                                designMat1, basisMat1,  lambda = NULL, nlam = 100, numCovs,
                                maxInt = 10^5,  epsilon = 1E-20, shrinkScale,accelrt=FALSE, truncation = TRUE){
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
  
  see = fitProxGradCppSeq(ulam, theta, stepSize, dat, basisMat0_tilda, n.k, Hp,
                          maxInt, epsilon, shrinkScale, accelrt, numCovs, designMat1_tilda, 
                          truncation,neg2loglikSat=start_fit$neg2loglik_sat)
  
  
  return(see)
  
  
}


