## Final function: DPMM 
## Input: data x, concentration parameter a, base distribution for alpha (mu1, sigma1)
##         for beta (mu2, sigma2), auxillary parameters number aux, iterarion times iter
## Output: A list as mentioned before

dpmm <- function(data, a = 20, mu1 = 0, sigma1 = 1, mu2 = 0, sigma2 = 1, aux = 10, iter = 10){
  source("Initialization.R")
  source("par_update.R")
  source("Rearrange.R")
  cl = Initialization(data, a, aux, mu1,sigma1,mu2,sigma2)
  for (t in 1:iter) {
    cl = Rearrange(cl, a, aux, mu1, sigma1, mu2, sigma2, data)
    cl = par_update(cl, data)
  }
  return(cl)
}


# test
# setwd("~/Drug_Discovery/BOISE/")
# library(rstan)
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
# rm(list = ls())
# load("clustering.RData")
# u <- foo$scaled.x
# dat <- 1*(u > .5 )
# cl = dpmm(dat,a = 20, aux = 10, iter = 2)

# s1 = apply(dat[which(cl$C == 4), ], 2, sum)
# sort(s1[which(s1>0)],decreasing = T)
# order(s1[which(s1>0)], decreasing = T)
# 
# s2 = apply(dat[which(cl$C == 5), ], 2, sum)
# sort(s2[which(s2>0)],decreasing = T)
# order(s2[which(s2>0)], decreasing = T)
# 
# s3 = apply(dat[which(cl$C == 6), ], 2, sum)
# sort(s3[which(s3>0)],decreasing = T)
# order(s3[which(s3>0)], decreasing = T)
# 
# rotate <- function(x) t(apply(x, 2, rev))
# image(rotate(dat[which(cl$C == 1),]))
# image(rotate(dat[which(cl$C == 2),]))
# image(rotate(dat[which(cl$C == 3),]))
