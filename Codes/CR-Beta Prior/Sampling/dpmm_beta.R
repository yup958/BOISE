## MCMC sampling for conjugate prior (Beta-binomial) clustering assignments
## Input: prior vector a,b, data matrix x0, warm-up step length warm, sample size iter, 
## gap between two samples step, divergence alpha 
## Output: A list cl_sample, including clustering assignments of size "iter * Start"

dpmm_beta <- function(a, b, x0, warm = 200, size = 10, iter = 10, step = 5, alpha = 2){
  source("Initial_beta.R")
  source("Update_beta.R")
  
  ## MCMC for random clustering
  n = dim(x0)[1]
  m = dim(x0)[2]
  KK = rep(0,iter)
  NN = matrix(0, iter, 2*n)
  CC = matrix(0,iter, n)
  # xx = matrix(0, iter, m)
  XX = rep(0, iter*size*m)
  dim(XX) = c(iter, size, m)
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

    ## sample for xi
    p = rep(0, cl$K + 1)
    p[cl$K + 1] = alpha / (n + alpha)
    p[1:cl$K] = sapply(1:cl$K,function(x){return(cl$N[x]/(n+alpha))})
    XX[i, , ] = t(as.matrix(sapply(1:size, function(s){
      classi = which(rmultinom(1, 1, p) == 1)
      if(classi == cl$K+1){
        posta = a
        postb = b
      } else{
        targets = which(cl$C == classi)
        posta = a
        postb = b
        for (t in targets) {
          posta = posta + x0[t, ]
          postb = postb + 1 - x0[t,]
        }
      }
      p = sapply(1:m, function(x){
        return(rbeta(1,posta[x],postb[x]))
      })
      new_xi = sapply(1:m, function(x){
        return(as.numeric(rbinom(1, 1, p = p[x])))
      })
      return(new_xi)
    })))
  }
  
  cl_sample = list(KK = KK, NN = NN, CC = CC, XX = XX)
  return(cl_sample)
}
# # Test
# system.time({
#   cl_sample = dpmm_beta(a,b,x0 = dat, warm = 500,size = 10,iter = 100, step = 10, alpha = 15)
# })
