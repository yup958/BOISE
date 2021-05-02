dpmm_beta <-
function(x0, alpha, beta, m0, burn_in = 200, sample_size = 10, thinning = 5){
#  source("Initial_beta.R")
#  source("Update_beta.R")
  
  ## Input check for x0, alpha/beta and m0
  m = nrow(x0)
  n = ncol(x0)
  for (i in 1:m) {
    for (j in 1:n) {
      if(!x0[i, j] %in% c(0,1)) stop("Invalid x0. x0 Should be binary matrix")
    }
  }
  if(length(alpha) != n | length(beta) != n)  stop("Length of prior parameter should be consistent")
  if(m0 <= 0) stop("Prior mass m0 should be greater than 0.")
  ## MCMC for random clustering
  
  KK = rep(0,sample_size)
  NN = matrix(0, sample_size, 2*m)
  CC = matrix(0,sample_size, m)
  ### Initialization
  cl = Initial_beta(x0, m0)
  ### burn-in stage
  for (i in 1:burn_in) {
    cl = Update_beta(cl, x0, alpha, beta, m0)
  }
    
  ### Sample clustering assignments
  for (i in 1:sample_size) {
    tmp = 0
    while (tmp < thinning) {
      cl = Update_beta(cl, x0, alpha, beta, m0)
      tmp = tmp + 1
    }
    KK[i] = cl$K
    NN[i, ] = cl$N
    CC[i, ] = cl$C
  }
  
  cl_sample = list(KK = KK, NN = NN, CC = CC)
  return(cl_sample)
}
