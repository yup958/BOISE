Boise_Aug <-
function(cl_sample, sample_size, interm_size, nT, alpha, beta, x0, m0, inform, nAdd){
  #source("clust_sum.R")
  #source("npel1.R")
  if (!require('parallel')) {
    stop('The package parallel was not installed')
  }
  P = clust_sum(cl_sample,x0,sample_size, alpha, beta)
  m = nrow(x0)
  n = ncol(x0)
  cl = cl_sample
  cl$XX = rep(0, sample_size * interm_size * n)
  dim(cl$XX) = c(sample_size, interm_size, n)
  
  ## Sample for x_i* with sample size "size"
  for (j in 1:sample_size) {
    K = cl$KK[j]
    p = rep(0, K + 1)
    p[K + 1] = m0 / (m + m0)
    p[1:K] = P[[j]][,n+1] / (m + m0)
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
  candidate = (1:ncol(x0))[-inform]
  while (step <= nAdd) {
    pel = rep(0,length(candidate))
    pel = unlist(mclapply(candidate, function(x){
      return(pel1_beta(cl, P, sample_size, interm_size, A = c(inform,x), nT,alpha, beta, x0, m0))},
      mc.cores = detectCores()))
    tmp = candidate[order(pel)[1]]
    inform = c(inform, tmp)
    step = step + 1
    candidate = candidate[-which(candidate == tmp)]
  }
  return(inform)
}
