### Dirichlet process updating for one iteration (based on Neal's 2000 Paper, Algorithm 3)
## Input: prior vectors a, b, divergence alpha, data, old cl list of {K, N, C}
## Output a new cl list {K, N, C}

Update_beta <- function(cl, dat, a, b, alpha = 2){
  n = dim(dat)[1]
  m = dim(dat)[2]
  new_K = cl$K
  new_N = cl$N
  new_C = cl$C
  
# Update labels
  for (i in 1:n) {
    if(new_N[new_C[i]] == 1){
      p = rep(0, new_K)
      for (k in 1:new_K) {
        if(new_N[k] == 0){
          next
        } else if(k == new_C[i]){
          p[k] = exp(log(alpha) + sum(dat[i,] * log(a / (a + b)), na.rm = T) + 
                       sum((1 - dat[i,]) * log(b / (a + b)), na.rm = T))
        } else{
          subdat = dat[which(new_C == k),]
          if(is.vector(subdat)){
            missing = which(is.na(subdat))
            ak = a + subdat
            ak[missing] = a[missing]
            bk = b + 1 - subdat
            bk[missing] = b[missing]
          } else{
            ak = a + colSums(subdat, na.rm = T)
            bk = b + colSums(1-subdat, na.rm = T)
          }
          q = ak / (ak + bk)
          p[k] = exp(log(new_N[k]) + sum(dat[i, ] * log(q), na.rm = T) + 
                       sum((1 - dat[i,]) * log(1 - q), na.rm = T))
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
      p[new_K + 1] = exp(log(alpha) + sum(dat[i,] * log(a / (a + b)), na.rm = T) + 
                           sum((1 - dat[i,]) * log(b / (a + b)), na.rm = T))
      for (k in 1:new_K) {
        if(new_N[k] == 0){
          next
        } else{
          subdat = dat[which(new_C == k),]
          if(is.vector(subdat)){
            missing = which(is.na(subdat))
            ak = a + subdat
            ak[missing] = a[missing]
            bk = b + 1 - subdat
            bk[missing] = b[missing]
          } else{
            ak = a + colSums(subdat,na.rm = T)
            bk = b + colSums(1-subdat,na.rm = T)
          }
          q = ak / (ak + bk)
          p[k] = exp(log(new_N[k]) + sum(dat[i, ] * log(q), na.rm = T) + 
                       sum((1 - dat[i,]) * log(1 - q), na.rm = T))
        }
      }
      p = p / sum(p)
      class = which(rmultinom(1,1,p) == 1)
      new_C[i] = class
      new_N[class] = new_N[class] + 1
      if(class == new_K + 1){new_K = new_K + 1}
    }
  }
  
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