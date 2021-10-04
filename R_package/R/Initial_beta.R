Initial_beta <-
function(x0, m0 = 2){
  m = nrow(x0)
  n = ncol(x0)
  N = rep(0, 2 * m)
  C = rep(0, m)
  K = 1
  
  N[1] = 1
  C[1] = 1
  for (i in 2:m) {
    p = rep(0, K + 1)
    p[K + 1] = m0 / (i - 1 + m0)
    p[1:K] = sapply(1:K,function(x){return(N[x]/(i-1+m0))})
    class = which(rmultinom(1, 1, p) == 1)
    C[i] = class
    N[class] = N[class] + 1
    if(class == K + 1){K = K + 1}
  }
  cl = list("K" = K, "N" = N, "C" = C)
  return(cl)
}
