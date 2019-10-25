## Input: Informer set A, Data matrix x0, Intermediate data xA, Informer size nA, Top size nT, 
## clustering assignment cl, prior a,b, divergence alpha
## Output: Top set T, PEL2 loss for x0, cl and xA, pel2, and posterior probability of P(xA|C,x0)

pel2_beta <- function (cl, x0, xA, nA, A, nT = 10, a, b,alpha){
  K = cl$K
  N = cl$N
  C = cl$C
  m = dim(x0)[2]
  p = rep(0, K + 1) #posterior probabilities of p(xA|ck,C,x0)
  theta = matrix(0, K + 1, m) # posterior expectations E(theta|x0,xA,ck,C)
  
  p[K + 1] = exp(sum(xA * log(a[A] / (a[A] + b[A]))) + sum((1 - xA) * log(b[A] / (a[A] + b[A]))))
  for (j in 1:m) {
    if(j %in% A){
      theta[K + 1, j] = (a[j] + xA[which(A == j)]) / (a[j] + b[j] + 1)
    } else{
      theta[K + 1, j] = a[j] / (a[j] + b[j])
    }
  }
  
  for (k in 1:K) {
    targets = which(C == k)
    a0 = a
    b0 = b
    for (t in targets) {
      a0 = a0 + x0[t, ]
      b0 = b0 + 1 - x0[t, ]
    }
    
    p[k] = exp(sum(xA * log(a0[A] / (a0[A] + b0[A]))) + sum((1 - xA) * log(b0[A] / (a0[A] + b0[A]))))
    for (j in 1:m) {
      if(j %in% A){
        theta[k, j] = (a0[j] + xA[which(A == j)]) / (N[k] + a0[j] + b0[j] + 1)
      } else{
        theta[k, j] = a0[j] / (N[k] + a0[j] + b0[j])
      }
    }
  }
  
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
  postls <- list(pel2 = 100, top = rep(0,nT), post_prob = post_prob)
  postls$top = order(post_theta, decreasing = T)[1:nT]
  postls$pel2 = sum(1 - post_theta[postls$top])
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
