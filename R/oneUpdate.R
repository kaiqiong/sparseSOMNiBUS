# One update for theta using backtracing line search
oneUpdate <- function(theta, stepSize, lambda1, dat, basisMat0, basisMat1, n.k, Hp, maxInt_lineSearch = 1000, numCovs,
                      shrinkScale,designMat1, theta_m, iter, accelrt){
  binom_out <- binomObject(theta,basisMat0, basisMat1, dat, n.k,numCovs, designMat1)
  gBinomLoss <- binom_out$gNeg2loglik
  
  new_out <- thetaUpdate(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, basisMat1, dat, designMat1, theta_m,
                         iter, accelrt)
  
  linearSupport <-  binom_out$neg2loglik - stepSize * sum(gBinomLoss * new_out$G_t_theta) + stepSize/2*sum(new_out$G_t_theta * new_out$G_t_theta)
  
  shrinkCondi <- new_out$binom_out_new$neg2loglik >linearSupport
  i = 1
  while(shrinkCondi & i < maxInt_lineSearch ){
    stepSize=stepSize*shrinkScale
    new_out <- thetaUpdate(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, basisMat1, dat, designMat1, theta_m,
                           iter, accelrt)
    # print(stepSize)
    linearSupport <-  binom_out$neg2loglik - stepSize * sum(gBinomLoss * new_out$G_t_theta) + stepSize/2*sum(new_out$G_t_theta * new_out$G_t_theta)
    shrinkCondi <- new_out$binom_out_new$neg2loglik >linearSupport
    i <- i +1
  }
  
  return(list(G_t_theta = new_out$G_t_theta, 
              theta_l_proximal=new_out$theta_l_proximal, 
              binom_out_new=new_out$binom_out_new,
              stepSize=stepSize,
              theta_l_proximal_sep= new_out$binom_out_new$theta.sep))
  
}




# return -2l(theta) and derivative of (-2l(theta))


# one theta update from the prox(theta_old - t gradient(theta_old))
thetaUpdate <- function(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, basisMat1, dat, designMat1, 
                        theta_m=NULL, iter,accelrt){
  if(accelrt){
    thetaInterm <- theta + (iter - 2)/(iter +1) *(theta-theta_m)
    theta_l <-  thetaInterm - stepSize * gBinomLoss # based only on binomial -2loglik function (Loss)
  }else{
    theta_l <-  theta - stepSize * gBinomLoss # based only on binomial -2loglik function (Loss)
  }
  thetaL.sep <- getSeparateTheta(theta_l, n.k, numCovs)
  theta_l_p_sep <- vapply(seq_len(numCovs+1), function(i){
    if(i ==1){
      thetaL.sep[[i]]
    }else{
      proximalOperator(t=stepSize, lambda1=lambda1, Hp=Hp, u_p=thetaL.sep[[i]] )
    }
  },rep(1, n.k))
  theta_l_proximal <- as.vector(theta_l_p_sep)
  
  G_t_theta <- (theta - theta_l_proximal)/stepSize
  
  binom_out_new <- binomObject(theta_l_proximal,basisMat0, basisMat1, dat, n.k, numCovs, designMat1)
  
  return(list(binom_out_new=binom_out_new,G_t_theta = G_t_theta, 
              theta_l_proximal=theta_l_proximal ))
}


thetaUpdateAcclrt <- function(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, basisMat1, dat, designMat1, theta_m, iter){
  
  thetaInterm <- theta + (iter - 2)/(inter +1) *(theta-theta_m)
  theta_l <-  thetaInterm - stepSize * gBinomLoss # based only on binomial -2loglik function (Loss)
  thetaL.sep <- getSeparateTheta(theta_l, n.k, numCovs)
  theta_l_p_sep <- vapply(seq_len(numCovs+1), function(i){
    if(i ==1){
      thetaL.sep[[i]]
    }else{
      proximalOperator(t=stepSize, lambda1=lambda1, Hp=Hp, u_p=thetaL.sep[[i]] )
    }
  },rep(1, n.k))
  theta_l_proximal <- as.vector(theta_l_p_sep)
  
  G_t_theta <- (theta - theta_l_proximal)/stepSize
  
  binom_out_new <- binomObject(theta_l_proximal,basisMat0, basisMat1, dat, n.k, numCovs, designMat1)
  
  return(list(binom_out_new=binom_out_new,G_t_theta = G_t_theta, 
              theta_l_proximal=theta_l_proximal ))
}

proximalOperator <- function(t, lambda1, Hp, u_p){
  
  threshold <-  as.numeric(1-t*lambda1/sqrt( sum( (u_p %*% Hp) * u_p)))
  if(threshold >= 0){
    return(threshold * u_p)
  }else{
    return(rep(0, length(u_p)))
  }
}