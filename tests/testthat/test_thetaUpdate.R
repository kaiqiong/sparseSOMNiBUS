

context("testing the function thetaUpdate(): proximal update with backtracking line search")

path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "ref_dat.RDS", sep = "")

#dat = readRDS(path_ref_data)
dat <- readRDS("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/sparseSOMNiBUS/tests/testthat/data/ref_dat.RDS")

dat <- dat[,1:10]
numCovs <- ncol(dat)-4


n.k = 10
lambda1 = 0

lambda2 = 0.1
stepSize=1
theta <- initTheta <- rep(0, n.k*(ncol(dat)-3))
shrinkScale=0.5

maxInt = 100
maxInt_lineSearch = 10
epsilon = 1E-6
printDetail = TRUE

accelrt = FALSE

out0 <- smoothConstructExtract(n.k, dat$Position, constrains = T)
out1 <- smoothConstructExtract(n.k, dat$Position, constrains = F)

basisMat0 <- out0$basisMat
basisMat1 <- out1$basisMat
smoOmega1 <- out1$smoothOmegaCr
sparOmega <- sparseOmegaCr(out1$myh, n.k, out1$matF) # the same for both intercept and non_intercept # Call a RCpp function
numCovs <- ncol(dat) - 4
designMat1 <- extractDesignMat1(numCovs, basisMat1, dat)

out <- binomObject(theta, basisMat0, basisMat1, dat, n.k, numCovs, designMat1)

# naive way to calculate the gradient

dim(designMat1[[1]])

all(designMat1[[1]][5,]==basisMat1[5,] * dat$X1[5])
all(designMat1[[6]][5,]==basisMat1[5,] * dat$X6[5])

mydesignMat <- cbind(basisMat0, designMat1[[1]], designMat1[[2]], 
                     designMat1[[3]],designMat1[[4]], 
                     designMat1[[5]], designMat1[[6]])

myg <- NULL
for( i in seq_len(length(theta))){
  myg <- c(myg, sum((dat$Meth_Counts - dat$Total_Counts * 0.5) * mydesignMat[,i] ))
}



test_that("binomObject returns the current values of -2loglik and its gradident", {
  expect_true(out$neg2loglik == -2*log(0.5)*sum(dat$Total_Counts))
  expect_true(all(-2*myg== out$gNeg2loglik))
})
