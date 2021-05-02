Boise <-
function(cl_sample, sample_size, interm_size, nA, nT, alpha, beta, x0, m0){
  #source("clust_sum.R")
  #source("npel1.R")
  if (!require('parallel')) {
    install.packages("parallel")
    library(parallel)
  }
  
  ### Input check
  m = nrow(x0)
  n = ncol(x0)
  for (i in 1:m) {
    for (j in 1:n) {
      if(!x0[i, j] %in% c(0,1)) stop("Invalid x0. x0 Should be binary matrix")
    }
  }
  if(length(alpha) != n | length(beta) != n)  stop("Length of prior parameter should be consistent")
  if(m0 <= 0) stop("Prior mass m0 should be greater than 0.")
  
  
  P = clust_sum(cl_sample, x0, sample_size, alpha, beta)
  cl = cl_sample
  cl$XX = rep(0, sample_size * interm_size * n)
  dim(cl$XX) = c(sample_size, interm_size, n)
  
  ## Sample for x_i* with sample size "size"
  for (j in 1:sample_size) {
    K = cl$KK[j]
    p = rep(0, K + 1)
    p[K + 1] = m0 / (m + m0)
    p[1:K] = P[[j]][ ,n+1] / (m + m0)
    cl$XX[j, , ] = t(as.matrix(sapply(1:interm_size, function(s){
      classi = which(rmultinom(1, 1, p) == 1)
      if(classi == K+1){
        post_theta = alpha / (alpha + beta)
      } else{
        post_theta = P[[j]][classi,1:n]
      }
      new_xi = sapply(1:n, function(x){
        return(as.numeric(rbinom(1, 1, p = post_theta[x])))
      })
      return(new_xi)
    })))
  }
  ## BOISE selection based on pel1
  step = 1
  pel1 = unlist(mclapply(1:ncol(x0), function(x){
    return(pel1_beta(cl, P, sample_size, interm_size, A = x, nT, alpha, beta, x0, m0))},
    mc.cores = detectCores()))
  
  tmp = order(pel1)[1]
  inform = tmp
  candidate = order(pel1)
  while (step < nA) {
    step = step +1
    candidate = candidate[-which(candidate == tmp)]
    pel = rep(0,length(candidate))
    pel = unlist(mclapply(candidate, function(x){
      return(pel1_beta(cl, P, sample_size, interm_size, A = c(inform,x), nT, alpha, beta, x0, m0))},
      mc.cores = detectCores()))
    tmp = candidate[order(pel)[1]]
    inform = c(inform, tmp)
  }
  return(inform)
}
