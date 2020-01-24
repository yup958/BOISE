## Input: Informer set A, Data matrix x0, Intermediate data xA, Informer size nA, Top size nT, 
## clustering assignment cl, prior a,b, divergence alpha
## Output: posterior expectation E(theta|C,x0,xA), and posterior probability P(xA|C,x0)

pel2_beta <- function (cl, x0, xA, nA, A, nT = 10, a, b,alpha){
  K = cl$K
  N = cl$N
  C = cl$C
  m = dim(x0)[2]
  n = dim(x0)[1]
  p = rep(0, K + 1) #posterior probabilities of p(xA|ck,C,x0)
  theta = matrix(0, K + 1, m) # posterior expectations E(theta|x0,xA,ck,C)
  
  theta[K+1,] = sapply(1:m, function(j){
    if(j %in% A){
      return((a[j] + xA[which(A == j)]) / (a[j] + b[j] + 1))
    } else{
      return(a[j] / (a[j] + b[j]))
    }
  })
  
  p[K + 1] = exp(sum(xA * log(a[A] / (a[A] + b[A]))) +
                   sum((1 - xA) * log(b[A] / (a[A] + b[A]))))
  
  p[1:K] = sapply(1:K, function(k){
    targets = which(C == k)
    a0 = a
    b0 = b
    for (t in targets) {
      a0 = a0 + x0[t, ]
      b0 = b0 + 1 - x0[t, ]
    }
    return(exp(sum(xA * log(a0[A] / (a0[A] + b0[A]))) + 
                 sum((1 - xA) * log(b0[A] / (a0[A] + b0[A])))))
  })
  
  theta[1:K,] = t(as.matrix(sapply(1:K, function(k){
    targets = which(C == k)
    a0 = a
    b0 = b
    for (t in targets) {
      a0 = a0 + x0[t, ]
      b0 = b0 + 1 - x0[t, ]
    }
    tmp = sapply(1:m, function(j){
      if(j %in% A){
        return((a0[j] + xA[which(A == j)]) / (N[k] + a0[j] + b[j] + 1))
      } else{
        return(a0[j] / (N[k] + a0[j] + b0[j]))
      }
    })
    return(tmp)
  })))
  
  ## P(xA|C,x0)
  m = rep(0, K + 1)
  m[1:K] = N[1:K] / (sum(N[1:K]) + alpha)
  m[K + 1] = alpha / (sum(N[1:K]) + alpha)
  post_prob = sum(m * p)
    
  ##E(theta|C,x0,xA)
  w = rep(0, K + 1)
  w[1:K] = N[1:K]
  w[K + 1] = alpha
  w = w * p
  w = as.matrix(w / sum(w))
  post_theta = t(w) %*% theta
  postls <- list(post_theta = post_theta, post_prob = post_prob)
  return(postls)
}

# Test
# x0 = dat
# xA = 0
# nA = 1
# A = 4
# nT = 10
# postls = pel2_beta(cl, x0, xA, nA, A, nT = 10, a, b, alpha = 2)
# postls$pel2
