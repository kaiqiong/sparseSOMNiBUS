# One update for theta using backtracing line search
oneUpdate <- function(theta, stepSize, lambda1, dat, basisMat0, n.k, Hp, numCovs,
                      shrinkScale,designMat1, theta_m, iter, accelrt){
  binom_out <- binomObject(theta,basisMat0, dat, n.k,numCovs, designMat1)
  gBinomLoss <- binom_out$gNeg2loglik
  
  
  new_out <- thetaUpdate(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, dat, designMat1, theta_m,
                         iter, accelrt)
  
  if(accelrt){
    tempDiff <- new_out$theta_l_proximal-new_out$thetaInterm
    linearSupport <-  new_out$binom_outInterm$neg2loglik + sum(new_out$binom_outInterm$gNeg2loglik *  tempDiff) +
      sum(tempDiff * tempDiff)/(2*stepSize)
  }else{
    linearSupport <-  binom_out$neg2loglik - stepSize * sum(gBinomLoss * new_out$G_t_theta) + stepSize/2*sum(new_out$G_t_theta * new_out$G_t_theta)
  }
  
  shrinkCondi <- new_out$binom_out_new$neg2loglik >linearSupport
  i = 1
  while(shrinkCondi){
    stepSize=stepSize*shrinkScale
    new_out <- thetaUpdate(stepSize=stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, dat, designMat1, theta_m,
                           iter, accelrt)
    # print(stepSize)
   # linearSupport <-  binom_out$neg2loglik - stepSize * sum(gBinomLoss * new_out$G_t_theta) + stepSize/2*sum(new_out$G_t_theta * new_out$G_t_theta)
    if(accelrt){
      tempDiff <- new_out$theta_l_proximal-new_out$thetaInterm
      linearSupport <-  new_out$binom_outInterm$neg2loglik + sum(new_out$binom_outInterm$gNeg2loglik *  tempDiff) +
        sum(tempDiff * tempDiff)/(2*stepSize)
    }else{
      linearSupport <-  binom_out$neg2loglik - stepSize * sum(gBinomLoss * new_out$G_t_theta) + stepSize/2*sum(new_out$G_t_theta * new_out$G_t_theta)
    }
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

#' @title One proximal gradient update
#'
#' @description One proximal gradient update
#' @param stepSize step size
#' @param theta estimate from the last step
#' @param gBinomLoss gradient of the loss function evalulated at
#' \code{theta}
#' @param n.k number of knots for all covariates (including intercept);
#' curretnly, we assume the same n.k for all covariates
#' @param lengthUniqueDataID number of samples in the data
#' @param numCovs number of covariates
#' @param lambda1 penalization parameter for the L2 norm
#' @param Hp penalization matrix, sum of smooth penalty and sparse penalty
#' given the value of lambda2
#' @param basisMat0 basis matrix for beta0(t)
#' @param basisMat1 basis matrix for betap(t), p = 1, 2, P
#' @param dat
#' @param designMat1
#' @param theta_m estimate from the previous step of theta
#' @param iter number of the current iteration
#' @param accelrt whether an accelerated approach is used or not
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{var.cov.alpha} var of alpha
#' \item \code{var.alpha.0} var of alpha0
#' \item \code{var.alpha.sep} var of alpha_p, p = 1, 2, P
#' }
#' @author Kaiqiong Zhao
#' @noRd
# one theta update from the prox(theta_old - t gradient(theta_old))
thetaUpdate <- function(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp, basisMat0, dat, designMat1, 
                        theta_m=NULL, iter,accelrt, truncation=TRUE){
  if(accelrt){
    thetaInterm <- theta + (iter - 2)/(iter +1) *(theta-theta_m)
    binom_outInterm <- binomObject(thetaInterm,basisMat0, dat, n.k,numCovs, designMat1, truncation)
    theta_l <-  thetaInterm - stepSize * binom_outInterm$gNeg2loglik # based only on binomial -2loglik function (Loss)
  }else{
    theta_l <-  theta - stepSize * gBinomLoss # based only on binomial -2loglik function (Loss)
  }
  thetaL.sep <- getSeparateTheta(theta_l, n.k, numCovs)
  
  #theta_l_p_sepAll <-
  #vapply(seq(nrow(gridIndex)), function(iii){
  #  vapply(seq_len(numCovs+1), function(i){
  #    if(i ==1){
  #      thetaL.sep[[i]]
  #    }else{
  #      proximalOperator(t=stepSize, lambda1=lambda1[gridIndex[iii,1]], Hp=HpAll[,,gridIndex[iii,2]], u_p=thetaL.sep[[i]] )
  #    }
  #  },rep(1, n.k))
  #}, matrix(1, nrow=n.k, ncol = numCovs+1) )
  
  theta_l_p_sep <- vapply(seq_len(numCovs+1), function(i){
    if(i ==1){
      thetaL.sep[[i]]
    }else{
      proximalOperator(t=stepSize, lambda1=lambda1, u_p=thetaL.sep[[i]] )
    }
  },rep(1, n.k))
  
  
  theta_l_proximal <- as.vector(theta_l_p_sep)

  #theta_l_proximal <- as.vector(theta_l_p_sepAll[,,jj])
  
  G_t_theta <- (theta - theta_l_proximal)/stepSize
  
  binom_out_new <- binomObject(theta_l_proximal,basisMat0, dat, n.k, numCovs, designMat1, truncation)
  
  if(accelrt){
    return(list(binom_out_new=binom_out_new,G_t_theta = G_t_theta, 
                theta_l_proximal=theta_l_proximal, thetaInterm=thetaInterm, binom_outInterm=binom_outInterm))
  }else{
    return(list(binom_out_new=binom_out_new,G_t_theta = G_t_theta, 
                theta_l_proximal=theta_l_proximal ))
  }
  

}


proximalOperator <- function(t, lambda1,  u_p){
  
  threshold <-  as.numeric(1-t*lambda1/sqrt( sum( u_p  * u_p)))
  if(threshold >= 0){
    return(threshold * u_p)
  }else{
    return(rep(0, length(u_p)))
  }
}