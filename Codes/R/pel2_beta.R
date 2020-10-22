pel2_beta <-
function (P, x0, xA, nA, A, nT = 10, a, b,alpha){
  K = nrow(P)
  m = ncol(x0)
  n = nrow(x0)
  p = rep(0, K + 1) #posterior probabilities of p(xA|ck,C,x0)
  theta = matrix(0, K + 1, m) # posterior expectations E(theta|x0,xA,ck,C)
  
  theta[K+1,] = a/(a+b)
  theta[K+1,A] = (a[A]+xA)/(a[A]+b[A]+1)
  
  # p(xA|i in ck, C, x0)
  p[K + 1] = exp(sum(xA * log(a[A] / (a[A] + b[A]))) + 
                   sum((1 - xA) * log(b[A] / (a[A] + b[A]))))
  p[1:K] = apply(P[,1:m],1,function(q){
    return(exp(sum(xA * log(q[A])) + sum((1-xA) * log(1 - q[A]))))
  })
  
  theta[1:K,] = P[,1:m]
  if(length(A) == 1){
    theta[1:K, A] = (P[,A] * (P[,m+1] + rep(a[1]+b[1], K)) + xA) /
      (P[,m+1] + rep(a[1]+b[1], K) + 1)
  } else{
    theta[1:K,A] = t(as.matrix(apply(P[,c(A,m+1)], 1, function(q){
      l = length(A)
      counts = q[1:l] * (q[l+1]+a[1]+b[1])
      counts = counts + xA
      return(counts / (q[l+1]+a[1]+b[1]+1))
    })))
  }
  theta[1:K, A] = t(matrix(rep(xA, K), nrow = nA, ncol = K))
  
  # P(xA|C,x0)
  weight = c(P[,m+1], alpha)
  weight = weight / (sum(weight))
  post_prob = sum(weight * p)
    
  ##E(theta|C,x0,xA)
  w = c(P[, m + 1], alpha)
  w = w * p
  w = as.matrix(w / sum(w))
  post_theta = t(w) %*% theta
  postls <- list(post_theta = post_theta, post_prob = post_prob)
  return(postls)
}
