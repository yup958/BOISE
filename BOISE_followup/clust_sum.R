clust_sum <-
function(cl_sample, x0, sample_size, alpha, beta){
  m = nrow(x0)
  n = ncol(x0)
  P = list(rep(0,sample_size))
  for (j in 1:sample_size) {
    K = cl_sample$KK[j]
    N = cl_sample$NN[j,1:K]
    tmp_cl = matrix(0,K,(n+1))
    tmp_cl[,n+1] = N
    tmp_cl[,1:n] = t(sapply(1:K,function(k){
      target = which(cl_sample$CC[j,] == k)
      if(length(target)==1){
        missing = which(is.na(x0[target,]))
        ak = alpha + x0[target,]
        bk = beta + (1 - x0[target,])
        ak[missing] = alpha[missing]
        bk[missing] = beta[missing]
        return(ak / (ak + bk))
      }else{
        ak = alpha + colSums(x0[target,], na.rm = T)
        bk = beta + colSums(1 - x0[target,], na.rm = T)
        return(ak / (ak + bk))
      }
    }))
    P[[j]] = tmp_cl
  }
  return(P)
}
