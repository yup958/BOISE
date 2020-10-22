### BOISE informer set selection, corresponding to Algorithm 2 in the paper
### Pre-step: Sample DPMM samples with dpmm_beta.R to get cl_sample object
### Input: list cl_sample, corresponding parameters a,b,iter,size,alpha; nA, nT; x0
### Output: Informer set A with size nA (via advanced adaptive selection)

Boise <- function(cl_sample, iter, size, nA, nT, a, b, x0, alpha){
  #source("clust_sum.R")
  #source("npel1.R")
  if (!require('parallel')) {
    stop('The package parallel was not installed')
  }
  P = clust_sum(cl_sample,x0,iter, a, b)
  n = nrow(x0)
  m = ncol(x0)
  cl = cl_sample
  cl$XX = rep(0, iter*size*m)
  dim(cl$XX) = c(iter, size, m)
  
  ## Sample for x_i* with sample size "size"
  for (j in 1:iter) {
    K = cl$KK[j]
    p = rep(0, K + 1)
    p[K + 1] = alpha / (n + alpha)
    p[1:K] = P[[j]][,m+1] / (n+alpha)
    cl$XX[j, , ] = t(as.matrix(sapply(1:size, function(s){
      classi = which(rmultinom(1, 1, p) == 1)
      if(classi == K+1){
        post_theta = a / (a + b)
      } else{
        post_theta = P[[j]][classi,1:m]
      }
      new_xi = sapply(1:m, function(x){
        return(as.numeric(rbinom(1, 1, p = post_theta[x])))
      })
      return(new_xi)
    })))
  }
  ## BOISE selection based on pel1
  step = 1
  pel1 = unlist(mclapply(1:dim(x0)[2], function(x){
    return(pel1_beta(cl, P, iter, size, A = x, nA = step, nT,a,b,x0, alpha))},
    mc.cores = detectCores()))
  
  tmp = order(pel1)[1]
  inform = tmp
  candidate = order(pel1)
  while (step < nA) {
    step = step +1
    candidate = candidate[-which(candidate == tmp)]
    pel = rep(0,length(candidate))
    pel = unlist(mclapply(candidate, function(x){
      return(pel1_beta(cl, P, iter, size, A = c(inform,x), nA = step,nT,a,b,x0, alpha))},
      mc.cores = detectCores()))
    tmp = candidate[order(pel)[1]]
    inform = c(inform, tmp)
  }
  return(inform)
}
