


sparseSmoothFit <- function(dat, n.k, stepSize=0.1,lambda1=0.002, lambda2=0.1, maxInt = 500,
                            epsilon = 1E-6, printDetail = TRUE, initTheta, shrinkScale=0.5,
                            accelrt = TRUE){
  
  # Insert some checks later
  
  initOut = extractMats(dat=dat,n.k=n.k)
  
  out = fitProxGrad(theta=initTheta, stepSize=stepSize,lambda1=lambda1,
              dat=dat, basisMat0= initOut$basisMat0, n.k=n.k, 
              sparOmega= initOut$sparOmega, lambda2=lambda2, 
              smoOmega1=initOut$smoOmega1,
              maxInt = maxInt, epsilon = epsilon, printDetail = printDetail,shrinkScale=shrinkScale,
              accelrt=accelrt, numCovs=initOut$numCovs, designMat1= initOut$designMat1, basisMat1=initOut$basisMat1)
  return(out)
  
}

# n.k number of knots
# pos a vector of genomic positions (no needed to ordered); i.e. the x values that the spline is built on

smoothConstructExtract <- function(n.k, pos, constrains = FALSE){
  if(!constrains){
    #see1 = mgcv::smooth.construct(mgcv::s(pos, k = n.k,
    #                                      fx = F, bs = "cr"), data = data.frame(pos), knots = NULL)
    see1= mgcv::smoothCon(mgcv::s(pos, k = n.k,
                                         fx = F, bs = "cr"), data = data.frame(pos), knots = NULL)[[1]] # see1 and see11 gives the same x, and the same S (within tolerance 10^-8)
    myh <- see1$xp[-1]-see1$xp[-n.k]   # smooth.construct directly order the x values
    matF <- matrix(see1$F, byrow = T, nrow = n.k)
    smoothOmegaCr <- see1$S[[1]]*see1$S.scale # this corresponds to the S matrix from smooth.construct !! same but numerically more stable
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

# Fit a model with sparseness-smoothness regularization for a fixed values of lambda1 and lambda2
fitProxGrad <- function(theta, stepSize,lambda1, dat, basisMat0, n.k, Hp, maxInt = 1000,
                        epsilon = 1E-6, printDetail = FALSE, shrinkScale,accelrt, numCovs, designMat1, basisMat1){
  
  #gridIndex <- expand.grid(seq(length(lambda1)), seq(dim(HpAll)[3]))
  
  #HpAll <- vapply(seq(length(lambda2)), function(ii){
  #  sparOmega + lambda2[ii]*smoOmega1 
  #}, sparOmega)
 # Hp <- sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0
 # Hp=(1-lambda2)*sparOmega + lambda2*smoOmega1
  iter <- 1

  out <- oneUpdate(theta, stepSize, lambda1, dat, basisMat0, n.k, Hp,  numCovs, shrinkScale, 
                         designMat1,theta_m=theta, iter=iter, accelrt)
  
  penTerms <- twoPenalties(out$theta_l_proximal_sep, lambda1, numCovs, n.k)
 
  Est.points <- rbind(theta, out$theta_l_proximal)
  geneGradL2 <- c(sqrt(sum(out$G_t_theta^2)))
  lossVals <- c(binomLoss = out$binom_out_new$neg2loglik, penTerms = penTerms)
  stepSizeVec <- c(stepSize, out$stepSize)
  while(sqrt(sum((out$theta_l_proximal -  theta)^2)) > epsilon & iter < maxInt){
    iter <- iter + 1
    theta_m <- theta
    theta <- out$theta_l_proximal  
    out <- oneUpdate(theta, stepSize, lambda1, dat, basisMat0, n.k, Hp,  numCovs, shrinkScale, 
                     designMat1,theta_m=theta_m, iter=iter, accelrt)
    
    penTerms <- twoPenalties(out$theta_l_proximal_sep,  lambda1, numCovs, n.k)
    Est.points <- rbind(Est.points, out$theta_l_proximal)
    geneGradL2 <- c(geneGradL2, sqrt(sum(out$G_t_theta^2)))
    lossVals <- rbind(lossVals, c(binomLoss = out$binom_out_new$neg2loglik, penTerms = penTerms))
    stepSizeVec <- c(stepSizeVec, out$stepSize)
    if(printDetail){
      print(iter)
      print(c(binomLoss = out$binom_out_new$neg2loglik, penTerms = penTerms))
    }
    #
  }

  zeroCovs <- unlist(lapply(out$theta_l_proximal_sep[-1], function(x){all(x==0)}))
  
  # No need to export the whold basisMat0 and basisMat1, the rows corresponding to diferent pos will be fine
  uniPos <- unique(dat$Position)
  uniPosID <- match(uniPos, dat$Position)
  
  # t
  
 # thetaEstSep <- lapply(out$theta_l_proximal_sep, function(x){as.vector(Linv %*% x) })
 # thetaEst <- unlist(thetaEstSep)
  
  return(list( thetaEst = out$theta_l_proximal,
              Est.points=Est.points, 
              geneGradL2=geneGradL2, 
              lossVals=lossVals, 
              stepSizeVec = stepSizeVec,
              thetaEstSep = out$theta_l_proximal_sep,
              basisMat0=basisMat0[uniPosID,], 
              basisMat1=basisMat1[uniPosID,],
              designMat1=designMat1, 
              lossSum = apply(lossVals, 1, sum),
              nzeros = sum(!zeroCovs),
              nzeroCovs = colnames(dat)[-c(1:4)][which(!zeroCovs)],
              uniPos = uniPos,
              sparOmega=sparOmega,
              smoOmega1=smoOmega1,
              numCovs=numCovs,
              zeroCovs=zeroCovs,
              lambda1 = lambda1,
              lambda2 = lambda2,
              pi_ij = out$binom_out_new$pi_ij,
              gNeg2loglik=out$binom_out_new$gNeg2loglik))
  
}

extractMats <- function(dat, n.k){
  out0 <- smoothConstructExtract(n.k, dat$Position, constrains = T)
  out1 <- smoothConstructExtract(n.k, dat$Position, constrains = F)
  
  basisMat0 <- out0$basisMat
  #basisMat0[,-1] <- scale(basisMat0[,-1])
  basisMat1 <- out1$basisMat
  smoOmega1 <- out1$smoothOmegaCr * nrow(dat)^2
  sparOmega <- sparseOmegaCr(out1$myh, n.k, out1$matF) # the same for both intercept and non_intercept # Call a RCpp function
  
  numCovs <- ncol(dat) - 4
  designMat1 <- extractDesignMat1(numCovs, basisMat1, dat)
  return(list(basisMat0=basisMat0,designMat1=designMat1,
              smoOmega1=smoOmega1,sparOmega=sparOmega,
              numCovs=numCovs, basisMat1=basisMat1))
}


  