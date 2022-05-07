Update_beta <-
function(cl, x0, alpha, beta, m0 = 2){
  m = nrow(x0)
  n = ncol(x0)
  new_K = cl$K
  new_N = cl$N
  new_C = cl$C
  
# Update labels
  for (i in 1:m) {
    if(new_N[new_C[i]] == 1){
      p = rep(0, new_K)
      for (k in 1:new_K) {
        if(new_N[k] == 0){
          next
        } else if(k == new_C[i]){
          p[k] = log(m0) + sum(x0[i,] * log(alpha / (alpha + beta)), na.rm = T) + 
                       sum((1 - x0[i,]) * log(beta / (alpha + beta)), na.rm = T)
        } else{
          subdat = x0[which(new_C == k),]
          if(is.vector(subdat)){
            missing = which(is.na(subdat))
            ak = alpha + subdat
            ak[missing] = alpha[missing]
            bk = beta + 1 - subdat
            bk[missing] = beta[missing]
          } else{
            ak = alpha + colSums(subdat, na.rm = T)
            bk = beta + colSums(1-subdat, na.rm = T)
          }
          q = ak / (ak + bk)
          p[k] = log(new_N[k]) + sum(x0[i, ] * log(q), na.rm = T) + 
                       sum((1 - x0[i,]) * log(1 - q), na.rm = T)
        }
      }
      p = p - max(p)
      p = exp(p) / sum(exp(p))
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
      p[new_K + 1] = log(m0) + sum(x0[i,] * log(alpha / (alpha + beta)), na.rm = T) + 
                           sum((1 - x0[i,]) * log(beta / (alpha + beta)), na.rm = T)
      for (k in 1:new_K) {
        if(new_N[k] == 0){
          next
        } else{
          subdat = x0[which(new_C == k),]
          if(is.vector(subdat)){
            missing = which(is.na(subdat))
            ak = alpha + subdat
            ak[missing] = alpha[missing]
            bk = beta + 1 - subdat
            bk[missing] = beta[missing]
          } else{
            ak = alpha + colSums(subdat,na.rm = T)
            bk = beta + colSums(1-subdat,na.rm = T)
          }
          q = ak / (ak + bk)
          p[k] = log(new_N[k]) + sum(x0[i, ] * log(q), na.rm = T) + 
                       sum((1 - x0[i,]) * log(1 - q), na.rm = T)
        }
      }
      p = p - max(p)
      p = exp(p) / sum(exp(p))
      class = which(rmultinom(1,1,p) == 1)
      new_C[i] = class
      new_N[class] = new_N[class] + 1
      if(class == new_K + 1){new_K = new_K + 1}
    }
  }
  
  new_cl = list("K" = 0, "N" = rep(0, 2 * m), "C" = rep(0,m))
  new_cl$K = length(which(new_N > 0))
  rank = order(new_N, decreasing = T)
  for (i in 1:new_cl$K) {
    r = rank[i]
    new_cl$N[i] = new_N[r]
    new_cl$C[which(new_C == r)] = i
  }
  
  return(new_cl)
}
