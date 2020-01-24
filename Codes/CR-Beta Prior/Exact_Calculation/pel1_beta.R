## Input: A sample of clustering assignments cl_sample, sample size iter, informer set A and nA,
## top size nT, prior a and b, data x0, divergence alpha
## Output: PEL1 of A

pel1_beta <- function(cl_sample, iter, A, nA, nT = 10, a, b, x0, alpha = 2){
  source("pel2.R")
  ### Compute PEL1
  KK = cl_sample$KK
  NN = cl_sample$NN
  CC = cl_sample$CC
  
  if(nA == 1){
      p = rep(0, 2)
      pel2 = rep(0, 2)
      theta = matrix(0, 2, dim(x0)[2])
      for (i in 1:iter) {
        tmp_cl = list(K = KK[i], N = NN[i, ], C = CC[i,])
        postls = pel2_beta(tmp_cl, x0, xA = 0, nA, A, nT, a, b, alpha)
        p[1] = p[1] + postls$post_prob
        theta[1, ] = theta[1, ] + postls$post_theta
        postls = pel2_beta(tmp_cl, x0, xA = 1, nA, A, nT, a, b, alpha)
        p[2] = p[2] + postls$post_prob
        theta[2,] = theta[2,] + postls$post_theta
      }
      theta = theta / iter
      p = p / iter
      p = p / sum(p)
      pel2[1] = sum(sort(1-theta[1,])[1:nT])
      pel2[2] = sum(sort(1-theta[2,])[1:nT])
  } else{
      p = rep(0, 2^nA)
      pel2 = rep(0, 2^nA)
      theta = matrix(0,2^nA,dim(x0)[2])
      
      ### Find all possible xA's
      X = matrix(0, 2^nA, nA)
      for (i in 1:(2^nA-1)) {
        s = i
        class(s) <- "binmode"
        s <- as.character(s)
        n = nchar(s)
        for (j in (nA-n+1):nA) {
          X[i+1,j] = as.integer(substr(s,j-nA+n,j-nA+n))
        }
      }
      
      
      for (i in 1:iter) {
        tmp_cl = list(K = KK[i], N = NN[i,], C = CC[i,])
        for (j in 1:2^nA) {
          postls = pel2_beta(tmp_cl, x0, xA = X[j, ], nA, A, nT, a, b, alpha)
          p[j] = p[j] + postls$post_prob
          theta[j, ] = theta[j, ] + postls$post_theta
        }
      }
      theta = theta / iter
      p = p / iter
      p = p / sum(p)
      pel2 = sapply(1:2^nA, function(x){
        return(sum(sort(1-theta[x,])[1:nT]))
      })
  }
  pel1 = sum(p * pel2)
  return(pel1)  
}

# Test
# x0 = dat
# A = 2
# nA = 1
# nT = 10
# a = apply(dat,2,sum)/100
# b = 2 - a
# iter = 10
# step = 5
# alpha = 2
# pel1 = pel1_beta(cl_sample, iter = 10, A,nA,nT,a,b,x0,alpha)
