
context("testing the RCpp function sparseOmegaCr.cpp")

#setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/snpSOMNiBUS/tests/testthat")

path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "ref_positions.RDS", sep = "")
#load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")
#saveRDS(pos, path_ref_data  )

pos = readRDS(path_ref_data)
n.k = 10
see1 = mgcv::smooth.construct(mgcv::s(pos, k = n.k,
                                      fx = F, bs = "cr"), data = data.frame(pos), knots = NULL)
myh <- see1$xp[-1]-see1$xp[-n.k] 
matF <- matrix(see1$F, byrow = T, nrow = n.k)
sparseOmegaCrR <- function(myh, K, matF){
  A11 <- A12 <- A22 <- matrix(0, nrow = K, ncol = K)
  for( i in 1:K){
    if(i==1){
      A11[i,i] <- myh[i]/3
      A12[i,i] <- -myh[i]^3/45
      A22[i,i] <- 4/315*myh[i]^5
    }
    if(i==K){
      A11[i,i] <- myh[i-1]/3
      A12[i,i] <- -myh[i-1]^3/45
      A22[i,i] <- 4/315*myh[i-1]^5
    }
    if(i >1 & i <K){
      A11[i,i] <- myh[i]/3 + myh[i-1]/3
      A12[i,i] <- -myh[i]^3/45-myh[i-1]^3/45
      A22[i,i] <- 4/315*myh[i]^5 + 4/315*myh[i-1]^5
    }
    if(i>1){
      A11[i, i-1] <- A11[i-1, i] <- myh[i-1]/6
      A12[i, i-1] <- A12[i-1, i] <- -7/360*myh[i-1]^3
      A22[i, i-1] <- A22[i-1, i] <- 31/15120*myh[i-1]^5
    }
  }
  return(A11 + t(matF)%*% t(A12) + A12 %*% matF + t(matF) %*% A22%*% matF)
}

cpp_out = sparseOmegaCr(myh, n.k, matF)
r_out = sparseOmegaCrR(myh, n.k, matF)
test_that("Rcpp returns the save sparse Omega for cubic regression as the naively coded R function sparseOmegaCrR ", {
  expect_true(is.matrix(cpp_out))
  expect_equal(nrow(cpp_out), n.k)
  expect_equal(ncol(cpp_out), n.k)
  expect_true(all(cpp_out - r_out <0.0000000001))
})
                              