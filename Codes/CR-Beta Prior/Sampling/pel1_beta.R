## Input: A sample of clustering assignments cl_sample, sample size iter, informer set A and nA,
## top size nT, prior a and b, data x0, divergence alpha
## Output: PEL1 of A

pel1_beta <- function(cl_sample, iter, size, A, nA, nT = 10, a, b, x0, alpha = 2){
  source("pel2.R")
  ### Compute PEL1
  KK = cl_sample$KK
  NN = cl_sample$NN
  CC = cl_sample$CC
  XX = cl_sample$XX
  
  pel1 = rep(0,iter)
  for (i in 1:iter) {
    tmp_cl = list(K = KK[i], N = NN[i,], C = CC[i,])
    XA = as.matrix(XX[i, ,A])
    bin = apply(XA,1,function(x){return(sum(x * 2^(0:(length(x)-1))))})
    XA = as.matrix(XA[order(bin),])
    bin = sort(bin)
    l = length(unique(bin))
    YA = matrix(0, l, nA)
    repli = rep(0, l)
    for (k in 1:l){
      tmp = which(bin == unique(bin)[k])
      YA[k,] = XA[tmp[1],]
      repli[k] = length(tmp)
    }
    ### Computed PEL1
    tmp_pel1 = sapply(1:l, function(x){
      post_theta = pel2_beta(tmp_cl, x0, xA = YA[x,], nA,A,nT,a,b,alpha)
      return(repli[x]*sum(sort(1-post_theta)[1:nT]))
      })
    pel1[i] = sum(tmp_pel1) / size
  }
  return(mean(pel1)) 
}

# # Test
# x0 = dat
# A = 2
# nA = 1
# nT = 10
# #a = apply(dat,2,sum)/100
# #b = 2 - a
# a = rep(mean(dat),366)
# b = 1 - a
# iter = 100
# step = 5
# alpha = 15
# system.time({
#   cl_sample = dpmm_beta(a,b,x0 = dat, warm = 100,iter = 100, step = 5, alpha = 15)
# })
# system.time({
#   pel1 = pel1_beta(cl_sample, iter = 100, size = 10,A = c(1,2),nA=2,nT,a,b,x0,alpha)
# })
