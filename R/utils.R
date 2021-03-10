

#'@param y meth_counts
#'@param x total_counts
#'@param designMat1  design matrix in the raw scale
#'@param basisMat0 design matrix in the raw scale
lambdaMax <- function(y, x, designMat1, Hp_inv, numCovs, basisMat0){
  
  nullFit <- glm(cbind(y, x-y)~-1+basisMat0, family = "binomial")
  
  mu <- nullFit$fitted.values
  
  diffvec = (y- x*mu)
  
  
 res=  vapply(1:numCovs, function(p){
    bp =  2* t(designMat1[[p]]) %*% diffvec
    sqrt(sum( (Hp_inv%*%bp) * bp))
    }, FUN.VALUE = c(0.0)
         )
  
  return(max(res))
  
}


getSatNeg2Loglik <- function(y,x, truncation = TRUE){
  eps <- 10 * .Machine$double.eps
  mu_s <- y/x
  if(truncation){
    mu_s[which(mu_s > 1 - eps)] <- 1 - eps
    mu_s[which(mu_s < eps)] <- eps
  }
  neg2loglik_sat <- ((-2)*sum(y*log(mu_s)+(x-y)*log(1-mu_s)))
  return(neg2loglik_sat)
}


#'@param ulam2 a sequence of lambda2
#'
getSeqLam1Hp <- function(lambda2,meth, total, lambda=NULL, nlam, sparOmega, smoOmega1, designMat1, basisMat0, hugeCont =100000 ){
  
  Hp <- (1-lambda2)*sparOmega + lambda2*smoOmega1 # H1 <- lambda2*smoOmega0
  
  #--- Matrix decomposition for Hp and calculate transformed design matrix
  
  L  = chol(Hp)
  Hinv = chol2inv(L)
  Linv = solve(L)
  #Linv = MASS::ginv(L)
  
  myp = (numCovs+1)*n.k
  lambda.min.ratio = ifelse(nrow(dat)<myp,0.01,0.0001)
  if (is.null(lambda)) {
    if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
    
    # compute lambda max: to add code here
    lambda_max <- lambdaMax(y=meth, x=total, 
                            designMat1=designMat1, Hp_inv=Hinv, numCovs=numCovs, basisMat0=basisMat0)
    
    
    # compute lambda sequence
    ulam <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                    length.out = nlam))
  } else { # user provided lambda values
    user_lambda = TRUE
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam = as.double(rev(sort(lambda)))
    nlam = as.integer(length(lambda))
  }
  
  ulam[1] = ulam[1]+hugeCont
  return(out = list(Linv=Linv, ulam = ulam ))
  
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


extractMats <- function(dat, n.k){
  out0 <- smoothConstructExtract(n.k, dat$Position, constrains = T)
  out1 <- smoothConstructExtract(n.k, dat$Position, constrains = F)
  
  basisMat0 <- out0$basisMat
  #basisMat0[,-1] <- scale(basisMat0[,-1])
  basisMat1 <- out1$basisMat
  smoOmega1 <- out1$smoothOmegaCr * nrow(dat)^2
  sparOmega <- sparseOmegaCr( out1$myh, n.k, out1$matF) # the same for both intercept and non_intercept # Call a RCpp function
  
  numCovs <- ncol(dat) - 4
  designMat1 <- extractDesignMat1(numCovs, basisMat1, dat)
  return(list(basisMat0=basisMat0,designMat1=designMat1,
              smoOmega1=smoOmega1,sparOmega=sparOmega,
              numCovs=numCovs, basisMat1=basisMat1))
}


# extract design matrix for beta1(t), ... betap(t)

extractDesignMat1 <- function(numCovs, basisMat1, dat ){
  
  lapply(seq_len(numCovs), function(i){(sweep(basisMat1, 1 ,dat[, paste0("X", i)], "*"  ))})
}
