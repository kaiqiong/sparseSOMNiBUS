# n.k number of knots
# pos a vector of genomic positions (no needed to ordered); i.e. the x values that the spline is built on

smoothConstructExtract <- function(n.k, pos, constrains = FALSE){
  if(!constrains){
    
  see1 = mgcv::smooth.construct(mgcv::s(pos, k = n.k,
                                        fx = F, bs = "cr"), data = data.frame(pos), knots = NULL)
  #see11= mgcv::smoothCon(mgcv::s(pos, k = n.k,
  #                                      fx = F, bs = "cr"), data = data.frame(pos), knots = NULL)[[1]] # see1 and see11 gives the same x, but different S, S.scale
  myh <- see1$xp[-1]-see1$xp[-n.k]   # smooth.construct directly order the x values
  matF <- matrix(see1$F, byrow = T, nrow = n.k)
  smoothOmegaCr <- see1$S[[1]]
  basisMat <- see1$X # nrow = length(pos), ncol = n.k, basis expansion matrix beta(pos) = basisMat %*% alpha
  }else
    {
  # the basis for beta0(t)  -- should be different from the one for beta1(t), because of additional constraints for beta0(t)
  see2 <- mgcv::smoothCon(mgcv::s(pos,k = n.k, fx = F, bs = "cr"),
                          data=data.frame(pos),knots=NULL, absorb.cons = TRUE)[[1]]
  basisMat <- cbind(rep(1, length(pos)), see2$X) 
  myh <- see2$xp[-1]-see2$xp[-n.k]   
  matF <- matrix(see2$F, byrow = T, nrow = n.k)
  smoothOmegaCr <- see2$S[[1]]
  }
  return(list(myh = myh, matF = matF, smoothOmegaCr = smoothOmegaCr, basisMat = basisMat))
}


# Assuming the effect of SNP all have the same n.k


# theta: basis coefficients for beta0(t),  beta1(t), beta2(t), ... betaP(t), a vector of length n.k *P -- Add this intercept or not?
# the first element of theta is an unpenalited constant c 
# length(theta) = n.k + n.k * numCovs

# numCovs: number of covariates 

# 
pos <- dat$Position
interceptNum <- 2
theta <- rnorm( (n.k-1) + n.k * n.snp, 3, 1 )
theta <- c(interceptNum, theta)
out0 <- smoothConstructExtract(n.k, pos, constrains = T)
out1 <- smoothConstructExtract(n.k, pos, constrains = F)

basisMat0 <- out0$basisMat
basisMat1 <- out1$basisMat

smoOmega0 <- out0$smoothOmegaCr
smoOmega1 <- out1$smoothOmegaCr

sparOmega <- sparseOmegaCr(out1$myh, n.k, out1$matF) # the same for both intercept and non_intercept

dat <- dat.use
numCovs <- n.snp

#sparseOmegaCr_out <- sparseOmegaCr(out0$myh, n.k-1, out0$matF)

getSeparateTheta <- function(theta, n.k, numCovs){
  theta.sep <- lapply(seq_len(numCovs+1), function(i){
    theta[ ((i-1)*n.k + 1): (i*n.k)]
  })
  return(theta.sep = theta.sep)
}

estimatePij <- function(theta,basisMat0, basisMat1, dat, n.k){
  numCovs <- ncol(dat) - 4
  
  if(length(theta)!= (n.k * (numCovs + 1))){
    stop("length(theta) is unequal to n.k * (numCovs + 1)")
  }

  theta.sep <- getSeparateTheta(theta, n.k, numCovs)
  # lp_p = beta(t) * Z =  (basisMatrix %*% alpha) *Z
  lp.sep <- vapply(seq_len(numCovs+1), function(i){ 
    if(i ==1){
      basisMat0 %*%  theta.sep[[i]]
    }else{
      (basisMat1 %*%  theta.sep[[i]]) *dat[, paste0("X", i-1)]
    }
  }, rep(1, nrow(dat))) # lp.sep = (beta0(t), beta1(t)*Z1, beta2(t) * Z2, ...)
  
  lp_ij <- apply(lp.sep, 1, sum)
  
  if(length(lp_ij) !=nrow(dat)){
    stop("length of lp_ij is unequal to nrow(dat)")
  }
  pi_ij <- 1/(1+exp(-lp_ij))
  return(list(pi_ij=pi_ij, theta.sep = theta.sep))

}

lambda2 = 2
lambda1 = 10


# return -2l(theta) and derivative of (-2l(theta))

binomObject <- function(theta,basisMat0, basisMat1, dat, n.k){ 
  
  numCovs <- ncol(dat) - 4
  # -2loglik for binomial
  estimatePijOut <- estimatePij(theta=theta, basisMat0 = basisMat0, basisMat1 = basisMat1, dat=dat, n.k = n.k)
  pi_ij <- estimatePijOut$pi_ij
  theta.sep <- estimatePijOut$theta.sep
  loglik_ij <- dat$Meth_Counts * log(pi_ij) +(dat$Total_Counts-dat$Meth_Counts)*log(1-pi_ij)
  
  neg2loglik <-  -2*sum(loglik_ij[is.finite(loglik_ij)])
  

  gPi_ij <- (dat$Meth_Counts - dat$Total_Counts*pi_ij) 
  # gradient for theta0
  see = sweep(basisMat0, 1, gPi_ij, "*" ) # first row * first element of gPi_ij ... 
  #all(see[1,] ==   basisMat0[1,]*gPi_ij[1])
  #all(see[2,] ==   basisMat0[2,]*gPi_ij[2] )
  g0_now = apply(see, 2, sum )
 
  # gradient for theta1, ... thetaP
  #basisMat1[1,] * dat[1, paste0("X", i)]
  #basisMat1[2,] * dat[2, paste0("X", i)]
  gRest <- vapply(seq_len(numCovs), function(i){
    designMatp <-  sweep(basisMat1, 1 ,dat[, paste0("X", i)], "*"  )
    apply(sweep(designMatp, 1, gPi_ij, "*"), 2, sum)
  }, rep(1, ncol(basisMat1)))
  
  gNeg2loglik <- -2*c(g0_now, as.vector(gRest))
  return(list(neg2loglik=neg2loglik, theta.sep = theta.sep,numCovs = numCovs,
              gNeg2loglik=gNeg2loglik))
}
twoPenalties <- function(theta.sep,smoOmega1, sparOmega, lambda1, lambda2, numCovs){
  
  #sparsePenalVector<- vapply(seq_len(numCovs+1), function(i){ 
  #  if(i ==1){
  #   (theta.sep[[i]][-1] %*%  sparOmega[2:n.k, 2:n.k]) %*% as.matrix(theta.sep[[i]][-1], ncol = 1)
     # see22 = sum(colSums((theta.sep[[i]] %*%  sparOmega) * theta.sep[[i]] ))
      
      
    #  autoplot(microbenchmark(see12 = (theta.sep[[i]] %*%  sparOmega) %*% as.matrix(theta.sep[[i]], ncol = 1),
    #                         see22 = sum(colSums((theta.sep[[i]] %*%  sparOmega) * theta.sep[[i]] )),
    #                          times = 50))
  #  }else{
  #    (theta.sep[[i]] %*%  sparOmega) %*% as.matrix(theta.sep[[i]], ncol = 1)
  #  }
  #}, 1) # lp.sep = (beta0(t), beta1(t)*Z1, beta2(t) * Z2, ...)
  
  #smoothPenalVector <- vapply(seq_len(numCovs+1), function(i){
  #  if(i==1){
  #    theta.sep[[i]][-1] %*%  smoOmega0 %*% as.matrix(theta.sep[[i]][-1], ncol = 1)
  #  }else{
  #    theta.sep[[i]] %*%  smoOmega1 %*% as.matrix(theta.sep[[i]], ncol = 1)
  #  }
  #}, 1) 
  
  Hp = sparOmega + lambda2*smoOmega1
  squaredIndPen <- vapply(seq_len(numCovs), function(i){
    sum(theta.sep[[i+1]] %*% Hp * theta.sep[[i+1]])
  },1)
  
  #all((squaredIndPen2-squaredIndPen)<0.000000001)
  #Hp = sparsePenalVector[-1] + lambda2 * smoothPenalVector[-1]   
  #penalTerm = lambda1 * sum(sparsePenalVector + lambda2 * smoothPenalVector)
  #squaredIndPen = sparsePenalVector[-1] + lambda2 * smoothPenalVector[-1] # t(theta) * Hp * theta
  penalTerm = lambda1 * sum(sqrt(squaredIndPen )) # no penalty on beta0(t)
  #sparsePenalVector[1] = 0
  #penalTerm = lambda1 * sum(sparsePenalVector + lambda2 * smoothPenalVector) # no sparse penalty on beta0(t) but with smooth penalty -- because the latter is readily available 
  return(list(penalTerm = penalTerm, Hp = Hp) )
}

twoPenalties <- function(theta.sep,Hp, lambda1, numCovs, n.k){
 # theta.sep <- getSeparateTheta(theta, n.k, numCovs)
  squaredIndPen <- vapply(seq_len(numCovs), function(i){
    sum(theta.sep[[i+1]] %*% Hp * theta.sep[[i+1]])
  },1)
  penalTerm = lambda1 * sum(sqrt(squaredIndPen )) # no penalty on beta0(t)
  
  return(penalTerm)
}
binom2PenalObject <- function(binom_out, penal_out){
  return(list(fullObj = binom_out$neg2loglik + penal_out$penalTerm))
}



binom_out <- binomObject(theta,basisMat0, basisMat1, dat, n.k)
penal_out <- twoPenalties(binom_out$theta.sep,smoOmega1, sparOmega, lambda1, lambda2, numCovs=binom_out$numCovs)
binom2PenalObject(binom_out, penal_out)
binom_out$neg2loglik
penal_out$penalTerm
#obj_out = binom2PenalObject(theta,basisMat0, basisMat1, dat,
#                  smoOmega0, smoOmega1, sparOmega, lambda1, lambda2, pos, n.k)


# the idea of the proximal operator (formula 4) is that one wants to minimize the penalty terms p(theta) at the same time we dont want to moving too much away from u
# prox_t(u) =  the proximal operator of the general group Lasso has an analytical solution

t=2
Hp = penal_out$Hp # depends on lambda1
u_p = rep(2, n.k)


proximalOperator <- function(t, lambda1, Hp, u_p){

 threshold <-  as.numeric(1-t*lambda1/sqrt( sum( (u_p %*% Hp) * u_p)))
 if(threshold >= 0){
   return(threshold * u_p)
 }else{
   return(rep(0, length(u_p)))
 }
}

proximalOperator(t=2, lambda1, Hp, u_p)

proximalOperator(t=20, lambda1, Hp, u_p)

#gradient for the -2loglik



# supplementary formula 4

theta.sep = binom_out$theta.sep
gNeg2loglik <- function(dat, pi_ij, numCovs, basisMat0, basisMat1){ 
  
  numCovs <-ncol(dat) - 4
  gPi_ij <- (dat$Meth_Counts - dat$Total_Counts*pi_ij) 
  
  # gradient for theta0
  see = sweep(basisMat0, 1, gPi_ij, "*" ) # first row * first element of gPi_ij ... 
      
  #all(see[1,] ==   basisMat0[1,]*gPi_ij[1])
  #all(see[2,] ==   basisMat0[2,]*gPi_ij[2] )
  g0_now = apply(see, 2, sum )
  #all((g0_now - g0)<0.0000001)
  
  #g0 <- 0
  #for( i in 1:nrow(basisMat0)){
  #  g0 = g0+ basisMat0[i,]*gPi_ij[i]
  #}
  
  # the design matrix should include the information of covariates
  
  # Step 1: obtain the design matrix basisMat1 * Z
  #basisMat1[1,] * dat[1, paste0("X", i)]
  #basisMat1[2,] * dat[2, paste0("X", i)]
  gRest <- vapply(seq_len(numCovs), function(i){
    
  designMatp <-  sweep(basisMat1, 1 ,dat[, paste0("X", i)], "*"  )
  apply(sweep(designMatp, 1, gPi_ij, "*"), 2, sum)
  
  }, rep(1, ncol(basisMat1)))
  
  return(-2*c(g0_now, as.vector(gRest)))
  
}
gradient_out <- gNeg2loglik(dat, binom_out$pi_ij, numCovs, basisMat0, basisMat1)





gdBinomLossUpdate <- function(theta, stepSize, gBinomLoss){
    return(theta - stepSize * gBinomLoss)
}



