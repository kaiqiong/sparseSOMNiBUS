# Assuming the effect of SNP all have the same n.k


# theta: basis coefficients for beta0(t),  beta1(t), beta2(t), ... betaP(t), a vector of length n.k *P -- Add this intercept or not?
# the first element of theta is an unpenalited constant c 
# length(theta) = n.k + n.k * numCovs

# numCovs: number of covariates 


# extract design matrix for beta1(t), ... betap(t)

extractDesignMat1 <- function(numCovs, basisMat1, dat ){
  
  lapply(seq_len(numCovs), function(i){(sweep(basisMat1, 1 ,dat[, paste0("X", i)], "*"  ))})
}


binomObject <- function(theta,basisMat0, dat, n.k,numCovs,designMat1, truncation = TRUE){ 
  
  # -2loglik for binomial
  estimatePijOut <- estimatePij(theta=theta, basisMat0 = basisMat0, designMat1 = designMat1, dat=dat, n.k = n.k, numCovs)
  pi_ij <- estimatePijOut$pi_ij
  eps <- 10 * .Machine$double.eps
  #if (any(pi_ij > 1 - eps) || any(pi_ij < eps)) 
  #  message("fitted probabilities numerically 0 or 1 occurred")
  #pi_ij_trc <- pi_ij
  if(truncation){
    pi_ij[which(pi_ij > 1 - eps)] <- 1 - eps
    pi_ij[which(pi_ij < eps)] <- eps
  }
  theta.sep <- estimatePijOut$theta.sep
  
  
  loglik_ij <- dat$Meth_Counts * log(pi_ij) +(dat$Total_Counts-dat$Meth_Counts)*log(1-pi_ij)
  
  #----------------------
  # if(truncation){
  #    loglik_ij[!is.finite(loglik_ij)] <- -10^20
  #  }
  #------------------------------
  #-2*sum(dbinom(x=dat$Meth_Counts, prob=pi_ij, size=dat$Total_Counts, log=T))
  neg2loglik <-  -2*sum(loglik_ij)
  
  
  
  # supplementary formula 4
  gPi_ij <- (dat$Meth_Counts - dat$Total_Counts*pi_ij) 
  # gradient for theta0
  see = sweep(basisMat0, 1, gPi_ij, "*" ) # first row * first element of gPi_ij ... 
  #all(see[1,] ==   basisMat0[1,]*gPi_ij[1]) # first observation
  #all(see[2,] ==   basisMat0[2,]*gPi_ij[2] )
  g0_now = apply(see, 2, sum )
  
  # gradient for theta1, ... thetaP
  #basisMat1[1,] * dat[1, paste0("X", i)]
  #basisMat1[2,] * dat[2, paste0("X", i)]
  gRest <- vapply(seq_len(numCovs), function(i){
    apply(sweep(designMat1[[i]], 1, gPi_ij, "*"), 2, sum)
  }, rep(1, ncol(designMat1[[1]])))
  
  gNeg2loglik <- -2*c(g0_now, as.vector(gRest))
  return(list(neg2loglik=neg2loglik, theta.sep = theta.sep,
              gNeg2loglik=gNeg2loglik, pi_ij=pi_ij))
}






getSeparateTheta <- function(theta, n.k, numCovs){
  theta.sep <- lapply(seq_len(numCovs+1), function(i){
    theta[ ((i-1)*n.k + 1): (i*n.k)]
  })
  return(theta.sep = theta.sep)
}

estimatePij <- function(theta,basisMat0, designMat1, dat, n.k, numCovs){
  
  if(length(theta)!= (n.k * (numCovs + 1))){
    message("length(theta) is unequal to n.k * (numCovs + 1)")
  }
  
  theta.sep <- getSeparateTheta(theta, n.k, numCovs)
  # lp_p = beta(t) * Z =  (basisMatrix %*% alpha) *Z
  lp.sep <- vapply(seq_len(numCovs+1), function(i){ 
    if(i ==1){
      basisMat0 %*%  theta.sep[[i]]
    }else{
      designMat1[[i-1]] %*%  theta.sep[[i]]
    }
  }, rep(1, nrow(dat))) # lp.sep = (beta0(t), beta1(t)*Z1, beta2(t) * Z2, ...)
  
  lp_ij <- apply(lp.sep, 1, sum)
  
  if(length(lp_ij) !=nrow(dat)){
    stop("length of lp_ij is unequal to nrow(dat)")
  }
  pi_ij <- 1/(1+exp(-lp_ij))
  return(list(pi_ij=pi_ij, theta.sep = theta.sep))
  
}


twoPenalties <- function(theta.sep,lambda1, numCovs, n.k){
  # theta.sep <- getSeparateTheta(theta, n.k, numCovs)
  squaredIndPen <- vapply(seq_len(numCovs), function(i){
    sum(theta.sep[[i+1]]* theta.sep[[i+1]])
  },1)
  penalTerm = lambda1 * sum(sqrt(squaredIndPen )) # no penalty on beta0(t)
  
  return(penalTerm)
}
#binom2PenalObject <- function(binom_out, penal_out){
#  return(list(fullObj = binom_out$neg2loglik + penal_out))
#}

# ssfit: an output from fitProxGradCpp()
optimcheckOld <- function(ssfit, lambda1, Hp, nk,eqDelta= 10^-5, uneqDelta=10^-5){
  
  thetaSep <- ssfit$thetaEstSep[-1]
  zeroCovs <- unlist(lapply( thetaSep, function(x){all(x==0)}))
  indOptimRes <- rep(NA, length(zeroCovs))
  
  zeroCovid <- which(zeroCovs)
  nonzeroCovid <- which(!zeroCovs)
 
  
  for(p in  nonzeroCovid ){
    
   ap = -ssfit$gNeg2loglik[(p*nk+1) :((p+1)*nk)]
    
   # ap =  2* t(designMat1[[p]]) %*% (dat$Meth_Counts-dat$Total_Counts*ssfit$pi_ij) 
    
  top = Hp %*% thetaSep[[p]]
  denom =sqrt(sum(top * thetaSep[[p]] )) 
  term2 = top/denom*lambda1
  
  indOptimRes[p]=sqrt(sum((ap - term2)^2)) <=eqDelta
    
  }
  
  for(p in zeroCovid){
    ap = -ssfit$gNeg2loglik[(p*nk+1) :((p+1)*nk)]
    
    temp1 = Hp %*% ap
   indOptimRes[p]= lambda1 -sqrt(sum(temp1*ap)) >= uneqDelta
  }
  return(indOptimRes)
  
  }
#Hp <- (1-lambda2)* sparOmega + lambda2*smoOmega1
#optimcheck(see2,lambda1=lambda1, Hp, nk = n.k)


#@thetaTildaSep  estimate of theta_tilda from fitProxGradCpp()  fit1$thetaEstSep
#@ gNeg2loglik  the gradient of negative twice log-likelihood, (w.r.t theta_tilda) fit1$gNeg2loglik

optimcheck <- function(thetaTildaSep, gNeg2loglik, lambda1, Hp, L, Linv, Hpinv, nk,eqDelta= 10^-5, uneqDelta=10^-5){
  
  # Note the output from fitProxGradCpp are theta_tilda not theta
  
  thetaSep <- thetaTildaSep[-1]
  
  thetaSep <- lapply(thetaSep, function(x){Linv %*% x })
  
  zeroCovs <- unlist(lapply( thetaSep, function(x){all(x==0)}))
  indOptimRes <- rep(NA, length(zeroCovs))
  
  zeroCovid <- which(zeroCovs)
  nonzeroCovid <- which(!zeroCovs)
  
  
  for(p in  nonzeroCovid ){
    
    bp = -gNeg2loglik[(p*nk+1) :((p+1)*nk)]; ap = t(L)%*% bp
    
    # ap =  2* t(designMat1[[p]]) %*% (dat$Meth_Counts-dat$Total_Counts*ssfit$pi_ij) 
    
    top = Hp %*% thetaSep[[p]]
    denom =sqrt(sum(top * thetaSep[[p]] )) 
    term2 = top/denom*lambda1
    
    indOptimRes[p]=sqrt(sum((ap - term2)^2)) <=eqDelta
    
  }
  
  for(p in zeroCovid){
    bp = -gNeg2loglik[(p*nk+1) :((p+1)*nk)]; ap = t(L)%*% bp
    
    #ap =  2* t(designMat1[[p]]) %*% (dat$Meth_Counts-dat$Total_Counts*ssfit$pi_ij) 
 
    temp1 = Hpinv %*% ap
    indOptimRes[p]= lambda1 -sqrt(sum(temp1*ap)) >= uneqDelta
  }
  return(indOptimRes)
  
}

