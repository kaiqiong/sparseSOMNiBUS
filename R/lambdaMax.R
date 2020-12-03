
lambdaMax <- function(y, x, designMat1, Hp, Hp_inv, numCovs, mu){
  
  
  diffvec = (y- x*mu)
 # b0 = 2* t(basisMat0) %*% diffvec
# sqrt(t(b0)%*% Hp %*% b0)
  
  res = NULL

  for(p in 1:numCovs){
  
  bp =  2* t(designMat1[[p]]) %*% diffvec 
  
  res = c(res, sqrt(sum( (Hp_inv%*%bp) * bp)))
  
  }
  
  return(max(res))
  
}
#lambdaMax(dat, designMat1, Hp, numCovs, mu)

#@param y meth_counts
#@param x total_counts


#Hp = (1-lambda2)*sparOmega*(nrow(dat)^2) + lambda2 * smoOmega1
