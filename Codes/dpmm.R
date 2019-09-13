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
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rm(list = ls())
load("clustering.RData")
u <- foo$scaled.x
dat <- 1*(u > .5 )
cl = dpmm(dat,a = 20, aux = 10, iter = 10)
