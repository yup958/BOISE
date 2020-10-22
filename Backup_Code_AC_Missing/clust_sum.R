### An intermediate function that summarizes the cl_sample to P,
##  A list of matrices, with P[[i]] a matrix contain all information of ith
##  cluster assignment. P[[i]] is K by (n+1) matrix. 
##  Each row for one cluster, last column as number of targets in the cluster

clust_sum <- function(cl_sample, dat, iter, a, b){
  n = nrow(dat)
  m = ncol(dat)
  P = list(rep(0,iter))
  for (j in 1:iter) {
    K = cl_sample$KK[j]
    N = cl_sample$NN[j,1:K]
    tmp_cl = matrix(0,K,(m+1))
    tmp_cl[,m+1] = N
    tmp_cl[,1:m] = t(sapply(1:K,function(k){
      target = which(cl_sample$CC[j,] == k)
      if(length(target)==1){
        missing = which(is.na(dat[target,]))
        ak = a + dat[target,]
        bk = b + (1 - dat[target,])
        ak[missing] = a[missing]
        bk[missing] = b[missing]
        return(ak / (ak + bk))
      }else{
        ak = a + colSums(dat[target,], na.rm = T)
        bk = b + colSums(1 - dat[target,], na.rm = T)
        return(ak / (ak + bk))
      }
    }))
    P[[j]] = tmp_cl
  }
  return(P)
}
