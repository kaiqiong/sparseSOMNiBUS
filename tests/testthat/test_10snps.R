# Samp.size: N  100
# Number of CpGs: 123 CpGs
# Number of SNPs: 20,60,100, 500, 1000, 5000
# Prop of significant SNPs:         5%
# Number:  1, 3, 5, 25, 50, 250

# Evaluate TPR




# Scenario 4:

# Simulate (Beta)-Binomial Outcomes

# the true methylated counts depends on 25 SNPs


# generate the non-zero betas

# 1. all with the exact same shape

# 2. varying shapes with either positive or negative effects on methylation




my.samp <- 100

n.snp <- 10


n.sig.snp <- round(n.snp *0.1)


library(magrittr)

BSMethSim_bbinom <-function(n, posit, theta.0, beta, phi, random.eff = F, mu.e=0,
                            sigma.ee=1,p0 = 0.003, p1 = 0.9, X, Z,binom.link="logit"){
  if( !is.matrix((Z)) ) message ("covariate Z is not a matrix")
  #  if( !is.matrix(beta) ) message ("the covariate effect parameter beta is not a matrix")
  
  if( !(nrow(X)==nrow(Z) & nrow(X) == n) ) message("Both X and Z should have n rows")
  if( !(ncol(X)==length(theta.0) & ncol(X) ==nrow (beta) & ncol(X)==length(posit) )) message ("The columns of X should be the same as length of beta theta.0 and posit; They all equals to the number of CpGs")
  if( ncol(beta)!= ncol(Z)) message("beta and Z should have the same dimentions")
  
  # the random effect term
  if(random.eff == T){
    my.e <- rnorm(n, mean=mu.e, sd = sqrt(sigma.ee))
  }else{
    my.e <- rep(mu.e, n)
  }
  
  my.theta <- t(sapply(1:n, function(i){
    theta.0 + rowSums(sapply(1:ncol(Z), function(j){Z[i,j] * beta[,j]})) + my.e[i]
  }))
  
  
  # Transform my.theta to my.pi for each (i, j)
  
  my.pi <- t(sapply(1:nrow(my.theta), function(i){
    #exp(my.theta[i,])/(1+exp(my.theta[i,]))
    binomial(link=binom.link)$linkinv(my.theta[i,])
  }))
  #  Generate S-ij based on the my.pi and my.Z
  #---------------------------------------------------#
  my.S <- my.pi
  for ( i in 1:nrow(my.S)){
    for ( j in 1:ncol(my.S)){
      #my.S[i,j] <- rbinom (1, size = X[i,j], prob= my.pi[i,j])
      my.S[i,j] <-VGAM::rbetabinom(1, size = X[i,j], prob = my.pi[i,j], rho = (phi[j]-1)/(X[i,j]-1) )
    }
  }
  #---------------------------------------------------#
  # Generate Y-ij based on the S-ij and the error rate (1-p1) and p0
  #---------------------------------------------------#
  my.Y <- my.S
  for ( i in 1:nrow(my.Y)){
    for ( j in 1:ncol(my.Y)){
      my.Y[i,j] <- sum(rbinom(my.S[i,j], size =1, prob=p1)) +
        sum(rbinom(X[i,j]-my.S[i,j], size = 1, prob=p0))
    }
  }
  out = list(S = my.S, Y = my.Y, theta = my.theta, pi = my.pi)
}

source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BSMethEM_April_update_scale_estp0p1_export_Q_known_error_with_phi_Mstep_quasi_correct_Fletcher.R")

#-------------------
# Step1 : Load a real data set 
#---------------------
load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")


# useful objects in this "BANK1data.RData" that I will use to simulate methylation data  

# methMat: matrix of methylated counts; rows are CpG sites, columns are samples
# totalMat: matrix of read-depth 
# pos: positions for the CpG sites in the methMat/totalMat
#my.pheno: the covariate files for the samles in methMat

#RAdat: this object is for running SOMNiBUS; it is not useful for simulation


#-----------------------------------------------------------------------------------------------
# Step2 : Specify the methylation patterns through functional parameters beta.0(t), beta.1(t), ...
#-------------------------------------------------------------------------------------------------

# I specify those values as output from real-data analysis using SOMNiBUS

load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1betas.RData")

# Load the functional parameters beta.0, beta.1, and beta.2 for intercept, effect of cell type and disease


beta.0 <- beta.0/4 + 2
#beta.0 <- beta.0-3.5



betas_non_zeros <- matrix(NA, nrow = length(beta.0), ncol = n.sig.snp)

for (i in 1:ncol(betas_non_zeros)){
  betas_non_zeros[,i] <- beta.1
}




beta_all <- cbind(beta.0, betas_non_zeros)


beta_all[,1] <- beta_all[,1]-30

#beta_all[,1] <- 0

# dispersion parameters are fixed as 1
phi <- rep(1,length(pos))


#-----------------------------------------------------------------------#


add_read_depth =0 # 



covs_use <- c("disease", "cell_type", "NullZ") 
#--------


#-------------------
# Step3 : Simulation
#---------------------

#---------------------------
# Step 3.1: simulate covariate matrix Z
#---------------------------

# the real-data covariate matrix is my.pheno; this step is to simulate Z based on my.pheno, with some simulation randomness



Z <- data.frame(matrix(NA, nrow= my.samp, ncol = n.snp))



set.seed(32314242)

mafs <- runif(n.snp,0.1, 0.5)
for ( i in 1:ncol(Z)){
  Z[,i] <- sample(c(0,1,2), size = my.samp, prob = c(mafs[i]^2, 2*mafs[i]*(1-mafs[i]), (1-mafs[i])^2), replace = T)
}


Z <-as.matrix(Z);rownames(Z)<- NULL

samp.Z <- Z


#---------------------------
# Step 3.2: Read depth matrix X
#---------------------------
# Build a read-depth matrix which sort of preserve the dependence structure in read-depths

my.X <- matrix(sample(0:1, my.samp*length(pos), replace = T), nrow = my.samp, ncol = length(pos))

ff = smooth.spline(pos, apply(totalMat, 1, median), nknots = 10)

spacial_shape <- round(predict(ff, pos)$y)
for ( i in 1:my.samp){
  my.X[i,] <- my.X[i,] + spacial_shape 
  #+ round(10*beta.0) + round(beta.1 * Z[i,1] + beta.2 * Z[i, 2])
} 

my.X <-  (my.X + add_read_depth)


#---------------------------
# Step 3.2: Simulate the methylated count matrix
#---------------------------

sim.dat<-BSMethSim_bbinom(n= my.samp, posit = pos, theta.0 =beta_all[,1], beta= beta_all[,-1], phi=phi, 
                          X = my.X, Z =Z[, 1:n.sig.snp],p0 = 0, p1 = 1,random.eff = F)



plot(pos, sim.dat$pi[1,], ylim = c(0,1))

for( i in 1:my.samp){
  points(pos, sim.dat$pi[i,], pch = 19, cex =0.5, col = i)
}



#-------------------------------------
# Generate loss function; proximal gradient descent; and etc.


#-- Step 1: generate the matrix of Omega1 and Omega2.


# n.k: K number of knots ---- 
# we choose the same k and same basis function B for all predictors 
# so, the same mat_omega2 for all p = 1, 2, ... P
length(pos)

# because curretnly, I didn't add smoothness penalty for the intercept, I will use n.k = 5 for the intercept
n.k = 20 # for the rest of Zs (non-intercept predictors)





#--- Organize the data before EM-smooth ---#
X <-my.X; Y <- sim.dat$Y 
samp.size <- nrow(Y); my.p <- ncol(Y)

dat.use <- data.frame(Meth_Counts=as.vector(t(Y)), 
                      Total_Counts=as.vector(t(X)), 
                      Position = rep(pos, samp.size),
                      ID = rep(1:samp.size, each=my.p))


covs_use <- colnames(samp.Z)
for( j in 1:length(covs_use)){
  dat.use <- data.frame(dat.use, rep(samp.Z[,j], each = my.p))
}
colnames(dat.use)[-c(1:4)] <- covs_use

dat.use <- dat.use[dat.use$Total_Counts>0,]

#my.span.dat<- data.frame(my.span.dat, null = sample(c(0,1), size = nrow(my.span.dat), replace = T))

#Z <- dat.use[,-c(1:4)]


# pre-set parameters
dat = dat.use
n.k = 10
numCovs = ncol(dat)-4
shrinkScale=1/2

n.k = 10
lambda1 = 0.002

lambda2 = 0.1
stepSize=3
theta <- initTheta <- rep(0, n.k*(ncol(dat)-3))
shrinkScale=0.5

maxInt = 100
maxInt_lineSearch = 10
epsilon = 1E-6
printDetail = TRUE

accelrt = FALSE

fit_out <- sparseSmoothFit(dat, n.k=n.k, stepSize=stepSize,lambda1=lambda1, lambda2=lambda2, maxInt = maxInt,
                           maxInt_lineSearch = 10, epsilon = 1E-6, printDetail = TRUE, initTheta=initTheta, shrinkScale = shrinkScale,
                           accelrt = accelrt)

saveRDS(fit_out, "/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/
        SOMNiBUS_SNP_selection/Rcpppackage/snpSOMNiBUS/tests/testthat/nSNP10noaccel.RDA")