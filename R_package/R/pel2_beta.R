pel2_beta <-
function (P, x0, xA, A, nT = 10, alpha, beta, m0){
  K = nrow(P)
  m = ncol(x0)
  n = nrow(x0)
  nA = length(A)
  lp = rep(0, K + 1) #posterior probabilities of log p(xA|ck,C,x0)
  theta = matrix(0, K + 1, m) # posterior expectations E(theta|x0,xA,ck,C)
  
  theta[K+1,] = alpha / (alpha + beta)
  theta[K+1,A] = (alpha[A]+xA)/(alpha[A]+beta[A]+1)
  
  # p(xA|i in ck, C, x0)
  lp[K + 1] = sum(xA * log(alpha[A] / (alpha[A] + beta[A]))) + 
                   sum((1 - xA) * log(beta[A] / (alpha[A] + beta[A])))
  if (K==1){
    lp[K] = sum(xA * log(P[A])) + sum((1-xA) * log(1 - P[A]))
  } else{
    lp[1:K] = apply(P[,1:m], 1, function(q){
      return(sum(xA * log(q[A])) + sum((1-xA) * log(1 - q[A])))
    })
  }
  
  theta[1:K,] = P[,1:m]
  if(length(A) == 1){
    theta[1:K, A] = (P[,A] * (P[,m+1] + rep(alpha[1]+beta[1], K)) + xA) /
      (P[,m+1] + rep(alpha[1]+beta[1], K) + 1)
  } else{
    theta[1:K,A] = t(as.matrix(apply(P[,c(A,m+1)], 1, function(q){
      l = length(A)
      counts = q[1:l] * (q[l+1]+alpha[1]+beta[1])
      counts = counts + xA
      return(counts / (q[l+1]+alpha[1]+beta[1]+1))
    })))
  }
  theta[1:K, A] = t(matrix(rep(xA, K), nrow = nA, ncol = K))
  
  # P(xA|C,x0)
  weight = c(P[,m+1], m0)
  log_weight = log(weight / (sum(weight)))
  log_weight = log_weight + lp
  c = max(log_weight)
  log_weight = log_weight - c
  log_post_prob = c + log(sum(exp(log_weight)))
    
  ##E(theta|C,x0,xA)
  lw = log(c(P[, m + 1], m0))
  lw = lw + lp
  lw = lw - max(lw)
  w = as.matrix(exp(lw) / sum(exp(lw)))
  post_theta = t(w) %*% theta
  postls <- list(post_theta = post_theta, log_post_prob = log_post_prob)
  return(postls)
}

