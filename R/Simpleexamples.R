

setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/src")
library(Rcpp)
sourceCpp("sparseOmegaCr.cpp")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/utils.R")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/sparseSmoothFit.R")
source("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/R/oneUpdate.R")

dat <- readRDS("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat/data/ref_dat.RDS")

dat <- dat[,1:10]

n.k = 10
lambda1 = 0

lambda2 = 0.1
stepSize=0.01 # Seems like I need a relative small step size, 

#stepSize = 1
theta <- initTheta <- rep(0, n.k*(ncol(dat)-3))
shrinkScale=0.5

maxInt = 5000
maxInt_lineSearch = 10
epsilon = 1E-6
printDetail = TRUE

accelrt = FALSE

fit_out <- sparseSmoothFit(dat, n.k=n.k, stepSize=stepSize,lambda1=lambda1, lambda2=lambda2, maxInt = maxInt,
                            maxInt_lineSearch = 10, epsilon = 1E-6, printDetail = TRUE, initTheta=initTheta, shrinkScale = shrinkScale,
                           accelrt = accelrt)
setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat")
saveRDS(fit_out, "nsnp6lambda0iter5000.RDA")
saveRDS(fit_out, "/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/
        SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat/nsnp6lambda0iter5000.RDA")

                "/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/
        SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat"

plot(fit_out$lossVals[,1]+fit_out$lossVals[,2], type = "l", 
     main= paste0("step size = ", stepSize, " MaxInt = ", maxInt, " Truncation = F"))

accelrt = T
fit_out_acc <- sparseSmoothFit(dat, n.k=n.k, stepSize=stepSize,lambda1=lambda1, lambda2=lambda2, maxInt = maxInt,
                           maxInt_lineSearch = 10, epsilon = 1E-6, printDetail = TRUE, initTheta=initTheta, shrinkScale = shrinkScale,
                           accelrt = accelrt)


plot(fit_out_acc$lossVals[,1]+fit_out_acc$lossVals[,2], type = "l")
