sparseSmoothFit <- function(dat, n.k, stepSize=0.1,lambda1=0.002, lambda2=0.1, maxInt = 500,
                            maxInt_lineSearch = 10, epsilon = 1E-6, printDetail = TRUE, initTheta, shrinkScale=0.5,
                            accelrt = TRUE){
  
  # Insert some checks later
  
  out0 <- smoothConstructExtract(n.k, dat$Position, constrains = T)
  out1 <- smoothConstructExtract(n.k, dat$Position, constrains = F)
  
  basisMat0 <- out0$basisMat
  basisMat1 <- out1$basisMat
  smoOmega1 <- out1$smoothOmegaCr
  sparOmega <- sparseOmegaCr(out1$myh, n.k, out1$matF) # the same for both intercept and non_intercept # Call a RCpp function
  
  numCovs <- ncol(dat) - 4
  designMat1 <- extractDesignMat1(numCovs, basisMat1, dat)
  
  out = fitProxGrad(theta=initTheta, stepSize=stepSize,lambda1=lambda1,
              dat=dat, basisMat0=basisMat0, basisMat1=basisMat1, n.k=n.k, sparOmega=sparOmega, lambda2=lambda2, smoOmega1=smoOmega1,
              maxInt = maxInt,
              maxInt_lineSearch = maxInt_lineSearch, epsilon = epsilon, printDetail = printDetail,shrinkScale=shrinkScale,
              accelrt=accelrt, numCovs, designMat1)
  return(out)
  
}

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


fitProxGrad <- function(theta, stepSize,lambda1, dat, basisMat0, basisMat1, n.k, sparOmega, lambda2, smoOmega1, maxInt = 1000,
                        maxInt_lineSearch = 1000, epsilon = 1E-6, printDetail = FALSE, shrinkScale,accelrt, numCovs, designMat1){
  
  
  Hp <- sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0
  iter <- 1

  out <- oneUpdate(theta, stepSize, lambda1, dat, basisMat0, basisMat1, n.k, Hp, maxInt_lineSearch, numCovs, shrinkScale, 
                         designMat1,theta_m=theta, iter=iter, accelrt)
  
  penTerms <- twoPenalties(out$theta_l_proximal_sep, Hp, lambda1, numCovs, n.k)
 
  Est.points <- rbind(theta, out$theta_l_proximal)
  geneGradL2 <- c(sqrt(sum(out$G_t_theta^2)))
  lossVals <- c(binomLoss = out$binom_out_new$neg2loglik, penTerms = penTerms)
  stepSizeVec <- c(stepSize, out$stepSize)
  while(sqrt(sum((out$theta_l_proximal -  theta)^2)) > epsilon & iter < maxInt){
    iter <- iter + 1
    theta_m <- theta
    theta <- out$theta_l_proximal  
    out <- oneUpdate(theta, stepSize, lambda1, dat, basisMat0, basisMat1, n.k, Hp, maxInt_lineSearch, numCovs, shrinkScale, 
                     designMat1,theta_m=theta_m, iter, accelrt)
    
    penTerms <- twoPenalties(out$theta_l_proximal_sep, Hp, lambda1, numCovs, n.k)
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
  return(list(thetaEst = out$theta_l_proximal,
              Est.points=Est.points, 
              geneGradL2=geneGradL2, 
              lossVals=lossVals, 
              stepSizeVec = stepSizeVec))
  
}

