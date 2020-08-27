
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "ref_dat.RDS", sep = "")
#load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")
#saveRDS(dat, path_ref_data  )

dat = readRDS(path_ref_data)

n.snp <- ncol(dat)-4

# pre-set parameters

n.k = 10
numCovs = ncol(dat)-4
shrinkScale=1/2


lambda2 = 2

lambda1 = 10

# step 1: Spline Basis Set up
# calculate matrices: sparseOmega, smoothOmega1, basisMat0 (intercept), basisMat1 (for rest of covariates)
# These matrices are fixed for fixed n.k and Position
out0 <- smoothConstructExtract(n.k, dat$Position, constrains = T)
out1 <- smoothConstructExtract(n.k, dat$Position, constrains = F)



basisMat0 <- out0$basisMat
basisMat1 <- out1$basisMat

smoOmega0 <- out0$smoothOmegaCr
smoOmega1 <- out1$smoothOmegaCr

setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/snpSOMNiBUS/src")
library(Rcpp)
sourceCpp("sparseOmegaCr.cpp")
sparOmega <- sparseOmegaCr(out1$myh, n.k, out1$matF) # the same for both intercept and non_intercept

Hp <- sparOmega + lambda2*smoOmega1

# Step 2: Given an initial value of: theta, lambda1, labmda2, update to get thetaNew
nterceptNum <- 2
theta <- rnorm( (n.k-1) + n.k * n.snp, 3, 1 )
theta <- c(interceptNum, theta)



maxInt_lineSearch <- maxInt <- 500
stepSize <- 1 # step size t


binom_out <- binomObject(theta,basisMat0, basisMat1, dat, n.k)
gBinomLoss <- binom_out$gNeg2loglik
#all(binom_out$gNeg2loglik==gradient_out)


numCovs <- ncol(dat) - 4
binom_out <- binomObject(theta,basisMat0, basisMat1, dat, n.k)
gBinomLoss <- binom_out$gNeg2loglik

binom_out$neg2loglik - stepSize * sum(out$binom_out$gNeg2loglik * out$G_t_theta) + t/2*sum(out$G_t_theta * out$G_t_theta)


thetaUpdate <- function(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp){
  theta_l <-  theta - stepSize * gBinomLoss # based only on binomial -2loglik function (Loss)
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
  
  binom_out_new <- binomObject(theta_l_proximal,basisMat0, basisMat1, dat, n.k)
  
  return(list(binom_out_new=binom_out_new,G_t_theta = G_t_theta, 
              theta_l_proximal=theta_l_proximal ))
}

oneUpdate <- function(theta, stepSize, lambda1, dat, basisMat0, basisMat1, n.k, Hp, maxInt_lineSearch = 1000){

  numCovs <- ncol(dat) - 4
  binom_out <- binomObject(theta,basisMat0, basisMat1, dat, n.k)
  gBinomLoss <- binom_out$gNeg2loglik
  
  new_out <- thetaUpdate(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp)
  
  linearSupport <-  binom_out$neg2loglik - stepSize * sum(gBinomLoss * new_out$G_t_theta) + t/2*sum(new_out$G_t_theta * new_out$G_t_theta)
  
  shrinkCondi <- binom_out_new$neg2loglik >linearSupport
  i = 1
  while(shrinkCondi & i < maxInt_lineSearch ){
    stepSize=stepSize*shrinkScale
    new_out <- thetaUpdate(stepSize, theta, gBinomLoss, n.k, numCovs, lambda1, Hp)
   # print(stepSize)
    linearSupport <-  binom_out$neg2loglik - stepSize * sum(gBinomLoss * new_out$G_t_theta) + t/2*sum(new_out$G_t_theta * new_out$G_t_theta)
    shrinkCondi <- new_out$binom_out_new$neg2loglik >linearSupport
    i <- i +1
  }
  
  return(list(G_t_theta = new_out$G_t_theta, 
              theta_l_proximal=new_out$theta_l_proximal, 
              binom_out_new=new_out$binom_out_new,
              stepSize=stepSize,
              theta_l_proximal_sep= new_out$binom_out_new$theta.sep))
 
}

lambda1 =0.001

fitProxGrad <- function(theta, stepSize,lambda1, dat, basisMat0, basisMat1, n.k, sparOmega, lambda2, smoOmega1, maxInt = 1000,
                        maxInt_lineSearch = 1000, epsilon = 1E-6, printDetail = FALSE){
  
  
  Hp <- sparOmega + lambda2*smoOmega1
  numCovs <- ncol(dat) - 4
  
  out <- oneUpdate(theta, stepSize, lambda1, dat, basisMat0, basisMat1, n.k, Hp, maxInt_lineSearch)
  penTerms <- twoPenalties(out$theta_l_proximal_sep, Hp, lambda1, numCovs, n.k)
  iter <- 1
  Est.points <- rbind(theta, out$theta_l_proximal)
  geneGradL2 <- c(sqrt(sum(out$G_t_theta^2)))
  lossVals <- c(binomLoss = out$binom_out_new$neg2loglik, penTerms = penTerms)
  stepSizeVec <- c(stepSize, out$stepSize)
  while(sqrt(sum((out$theta_l_proximal -  theta)^2)) > epsilon & iter < maxInt){
    iter <- iter + 1
    theta <- out$theta_l_proximal
    out <- oneUpdate(theta, stepSize, lambda1, dat, basisMat0, basisMat1, n.k, Hp, maxInt_lineSearch)
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


pgOut <- fitProxGrad(theta=rep(0, n.k*(n.snp+1)), stepSize=0.1,lambda1=0.002,
            dat, basisMat0, basisMat1, n.k, sparOmega, lambda2, smoOmega1,
            maxInt = 50,
            maxInt_lineSearch = 50, epsilon = 1E-6, printDetail = TRUE)

plot(pgOut[[4]][,1]+pgOut[[4]][,2])

plot(apply(lossVals, 1, sum), type = "l")


save.image(file = "nSNP1000.RData")



