### Dirichlet process updating for one iteration (based on Neal's 2000 Paper, Algorithm 3)
## Input: prior vectors a, b, divergence alpha, data, old cl list of {K, N, C}
## Output a new cl list {K, N, C}

Update_beta <- function(cl, dat, a, b, alpha = 2){
  n = dim(dat)[1]
  m = dim(dat)[2]
  new_K = cl$K
  new_N = cl$N
  new_C = cl$C
  # new_xi = cl$xi
  
  
  for (i in 1:n) {
    if(new_N[new_C[i]] == 1){
      p = rep(0, new_K)
      for (k in 1:new_K) {
        if(new_N[k] == 0){
          next
        } else if(k == new_C[i]){
          p[k] = exp(log(alpha) + sum(dat[i,] * log(a / (a + b))) + 
                       sum((1 - dat[i,]) * log(b / (a + b))))
        } else{
          targets = which(new_C == k)
          tmp = a
          for (t in targets) {
            tmp = tmp + dat[t, ]
          }
          q = tmp / (a + b + new_N[k])
          p[k] = exp(log(new_N[k]) + sum(dat[i, ] * log(q)) + 
                       sum((1 - dat[i,]) * log(1 - q)))
        }
      }
      p = p / sum(p)
      class = which(rmultinom(1,1,p) == 1)
      if(class != new_C[i]){
        new_K = new_K - 1
        new_N[new_C[i]] = 0
        new_N[class] = new_N[class] + 1
        new_C[i] = class
      }
    } else{
      new_N[new_C[i]] = new_N[new_C[i]] - 1
      new_C[i] = 0
      p = rep(0, new_K + 1)
      p[new_K + 1] = exp(log(alpha) + sum(dat[i,] * log(a / (a + b))) + 
                           sum((1 - dat[i,]) * log(b / (a + b))))
      for (k in 1:new_K) {
        if(new_N[k] == 0){
          next
        } else{
          targets = which(new_C == k)
          tmp = a
          for (t in targets) {
            tmp = tmp + dat[t, ]
          }
          q = tmp / (a + b + new_N[k])
          p[k] = exp(log(new_N[k]) + sum(dat[i, ] * log(q)) + 
                       sum((1 - dat[i,]) * log(1 - q)))
        }
      }
      p = p / sum(p)
      class = which(rmultinom(1,1,p) == 1)
      new_C[i] = class
      new_N[class] = new_N[class] + 1
      if(class == new_K + 1){new_K = new_K + 1}
    }
  }
  
  # ## sample for xi
  # p = rep(0, new_K + 1)
  # p[new_K + 1] = alpha / (n + alpha)
  # p[1:new_K] = sapply(1:new_K,function(x){return(new_N[x]/(n+alpha))})
  # classi = which(rmultinom(1, 1, p) == 1)
  # if(classi == new_K+1){
  #   posta = a
  #   postb = b
  # } else{
  #   targets = which(new_C == classi)
  #   posta = a
  #   postb = b
  #   for (t in targets) {
  #     posta = posta + dat[t, ] 
  #     postb = postb + 1 - dat[t,]
  #   }
  # }
  # p = sapply(1:m, function(x){
  #   return(rbeta(1,posta[x],postb[x]))
  # })
  # new_xi = sapply(1:m, function(x){
  #   return(as.numeric(rbernoulli(1, p = p[x])))
  # })
  
  # new_cl = list("K" = 0, "N" = rep(0, 2*n), "C" = rep(0,n), "xi" = new_xi)
  new_cl = list("K" = 0, "N" = rep(0, 2 * n), "C" = rep(0,n))
  new_cl$K = length(which(new_N > 0))
  rank = order(new_N, decreasing = T)
  for (i in 1:new_cl$K) {
    r = rank[i]
    new_cl$N[i] = new_N[r]
    new_cl$C[which(new_C == r)] = i
  }
  
  return(new_cl)
}


# #Test
# load("clustering.RData")
# load("pkis1.rda")
# u <- foo$scaled.x
# dat <- 1*(u > .5 )
# rm(u)
# rm(foo)
# dat <- t(apply(pkis1, 1, function(x){
#   thres = mean(x) + 2 * sd(x)
#   return(as.numeric(x>thres))
# }))
# #a = apply(dat,2,sum)/100
# #b = 2 - a
# a = rep(mean(dat),366)
# b = 1-a
# cl = Initial_beta(dat,a,b,alpha = 20)
# cl = Update_beta(cl, dat,a,b, alpha = 20)
# for (i in 1:200) {
#   cl = Update_beta(cl, dat,a,b, alpha = 20)
# }
# sum(cl$xi)
# which(cl$xi == 1)
