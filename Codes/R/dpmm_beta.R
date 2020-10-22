dpmm_beta <-
function(a, b, x0, warm = 200, iter = 10, step = 5, alpha = 2){
#  source("Initial_beta.R")
#  source("Update_beta.R")
  
  ## MCMC for random clustering
  n = dim(x0)[1]
  m = dim(x0)[2]
  KK = rep(0,iter)
  NN = matrix(0, iter, 2*n)
  CC = matrix(0,iter, n)
  ### Initialization
  cl = Initial_beta(x0, a, b, alpha)
  ### burn-in stage
  for (i in 1:warm) {
    cl = Update_beta(cl, dat = x0, a, b, alpha)
  }
    
  ### Sample clustering assignments
  for (i in 1:iter) {
    tmp = 0
    while (tmp < step) {
      cl = Update_beta(cl, dat = x0, a, b, alpha)
      tmp = tmp + 1
    }
    KK[i] = cl$K
    NN[i, ] = cl$N
    CC[i, ] = cl$C
  }
  
  cl_sample = list(KK = KK, NN = NN, CC = CC)
  return(cl_sample)
}
