### BOISE selection for augmentation of an existing informer set, usually used in HPC.
### Pre-step: Sample DPMM samples of cl_sample. 
###           A previous selected smaller informer set inform.
### Input: list cl_sample, corresponding parameters a,b,iter,size,alpha; nT; x0
###         Pre-selected informer set inform, number of informers to add nAdd
### Output: Informer set A with size n+nAdd (via advanced adaptive selection)

Boise_Aug <- function(cl_sample, iter, size, nT, a, b, x0, alpha, inform, nAdd){
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
  n = length(inform)
  candidate = (1:dim(x0)[2])[-inform]
  while (step <= nAdd) {
    pel = rep(0,length(candidate))
    pel = unlist(mclapply(candidate, function(x){
      return(pel1_beta(cl, P, iter, size, A = c(inform,x), nA = n+step,nT,a,b,x0, alpha))},
      mc.cores = detectCores()))
    tmp = candidate[order(pel)[1]]
    inform = c(inform, tmp)
    step = step + 1
    candidate = candidate[-which(candidate == tmp)]
  }
  return(inform)
}