Initial_beta <-
function(dat, a, b, alpha = 2){
  n = dim(dat)[1]
  m = dim(dat)[2]
  N = rep(0, 2 * n)
  C = rep(0, n)
  K = 1
  
  N[1] = 1
  C[1] = 1
  for (i in 2:n) {
    p = rep(0, K + 1)
    p[K + 1] = alpha / (i - 1 + alpha)
    p[1:K] = sapply(1:K,function(x){return(N[x]/(i-1+alpha))})
    class = which(rmultinom(1, 1, p) == 1)
    C[i] = class
    N[class] = N[class] + 1
    if(class == K + 1){K = K + 1}
  }
  cl = list("K" = K, "N" = N, "C" = C)
  return(cl)
}
