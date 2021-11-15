post_theta_j <- function(P, x0, xA, A, alpha, beta, m0){
  ### calculate E(mean_j theta_ij | x0, xA, C)
  K = nrow(P)
  m = ncol(x0)
  n = nrow(x0)
  nA = length(A)
  p = rep(0, K + 1) #posterior probabilities of p(xA|ck,C,x0)
  theta = matrix(0, K + 1, m) # posterior expectations E(theta|x0,xA,ck,C)
  
  ### p(xA|i in ck, C, x0)
  p[K + 1] = exp(sum(xA * log(alpha[A] / (alpha[A] + beta[A]))) + 
                   sum((1 - xA) * log(beta[A] / (alpha[A] + beta[A]))))
  p[1:K] = apply(P[,1:m],1,function(q){
    return(exp(sum(xA * log(q[A])) + sum((1-xA) * log(1 - q[A]))))
  })
  
  ### E(mean theta_ij | i* in ck, C, x0, xA)
  
  # if i* in new cluster:
  # sum of post expectations of theta_ij (or phi_kj) in each cluster
  post_exp_sum = P[ , 1:m] * P[ , (m+1)] 
  # Take the weighted average. Do not use colMeans!
  theta[K+1,] = colSums(post_exp_sum) / n
  
  # else:
  for(k in 1:K){
    post_exp_sum = P[ , 1:m]
    post_exp_sum[k, A] = (post_exp_sum[k, A] * (P[k, m+1] + 1) + xA) / (P[k, m+1] + 2)
    post_exp_sum = post_exp_sum * P[ , (m+1)]
    theta[k,] = colSums(post_exp_sum) / n
  }
  
  ##E(mean theta_j|C,x0,xA)
  w = c(P[, m + 1], m0)
  w = w * p
  w = as.matrix(w / sum(w))
  post_exp_mean = t(w) %*% theta
  return(post_exp_mean)
}
