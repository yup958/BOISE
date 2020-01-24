### Innitialization based on Chinese Restaurant Process for conjugate prior DPMM
## Input: number of targets to be clustered n, divergency parameter alpha
## Output: A list of {K, N, C,xi}. K = number of clusters, N = #of targets in each cluster, 
## C = cluster assignment (including new target), x_i* = sampled outcome on new target.

Initial_beta <- function(dat, a, b, alpha = 15){
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
  
  
  # ## sample for xi*
  # p = rep(0, K + 1)
  # p[K + 1] = alpha / (n + alpha)
  # p[1:K] = sapply(1:K,function(x){return(N[x]/(n+alpha))})
  # classi = which(rmultinom(1, 1, p) == 1)
  # if(classi == K + 1){
  #   posta = a
  #   postb = b
  # } else{
  #     targets = which(C == classi)
  #     posta = a
  #     postb = b
  #     if(length(targets) == 1){
  #       posta = posta + dat[targets,]
  #       postb = postb + 1 - dat[targets,]
  #     } else{
  #       posta = posta + apply(dat[targets, ], 2, sum)
  #       postb = postb + apply(1 - dat[targets, ], 2, sum)
  #     }
  # }
  # p = sapply(1:m, function(x){
  #   return(rbeta(1,posta[x],postb[x]))
  # })
  # xi = sapply(p, function(x){
  #   return(as.numeric(rbinom(1, 1, x)))
  # })
  # 
  # cl = list("K" = K, "N" = N, "C" = C, "xi" = xi)
  cl = list("K" = K, "N" = N, "C" = C)
  return(cl)
}

# 
#Test
# load("clustering.RData")
# u <- foo$scaled.x
# dat <- 1*(u > .5 )
# rm(u)
# rm(foo)
# a = rep(mean(dat),366)
# b = 1-a
# cl = Initial_beta(dat,a,b,10)
