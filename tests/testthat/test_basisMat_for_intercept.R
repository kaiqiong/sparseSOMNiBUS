
context("generating basis matrix under constrains")

# this test function can be removed later

#setwd("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_SNP_selection/Rcpppackage/snpSOMNiBUS/tests/testthat")

path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "ref_positions.RDS", sep = "")
#load("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1data.RData")
#saveRDS(pos, path_ref_data  )

pos = readRDS(path_ref_data)
n.k = 10
temp = sample(c(20:120), size = length(pos),replace = TRUE)

dat_check <- data.frame(pos = pos, X = temp, Y = temp-10 )

samp_gam <- mgcv::gam(cbind(Y, X-Y)~s(pos, k = n.k, fx = F, bs ="cr"), data = dat_check, family = binomial)


basisMat0 <- model.matrix(samp_gam)

sm <- mgcv::smoothCon(mgcv::s(pos,k = n.k, fx = F, bs = "cr"),data=data.frame(pos),knots=NULL, absorb.cons = TRUE)[[1]]


#----- test function smoothConstructOut()

const_out0 <- smoothConstructExtract(n.k, pos, constrains = TRUE)
const_out1 <- smoothConstructExtract(n.k, pos, constrains = FALSE)

test_that("smoothCon can generate basis under constrains", {
  expect_true(all(basisMat0 == cbind(rep(1, length(pos)), sm$X)))
  expect_equal(ncol(const_out1$basisMat), ncol(const_out0$basisMat))
  expect_equal(const_out0$basisMat[,1], rep(1, length(pos)))
  
})
                              