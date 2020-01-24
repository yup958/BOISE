### Informer set selection function version 1
### Pre-step: Sample DPMM samples of cl_sample
### Input: list cl_sample, corresponding parameters a,b,iter,size,alpha; nA, nT; x0
### Output: Informer set A with size nA (via advanced adaptive selection)

inform_beta1 <- function(cl_sample, iter, size, nA, nT, a, b, x0, alpha){
  source("pel1_beta.R")
  step = 1
  pel1 = sapply(1:dim(x0)[2], function(x){
    return(pel1_beta(cl_sample, iter, size, A = x, nA = step, nT,a,b,x0, alpha))})
  tmp = order(pel1)[1]
  inform = tmp
  candidate = order(pel1)
  while (step < nA) {
    step = step +1
    candidate = candidate[-which(candidate == tmp)]
    pel = rep(0,length(candidate))
    pel = sapply(candidate, function(x){
      return(pel1_beta(cl_sample, iter, size, A = c(inform,x), nA = step,nT,a,b,x0, alpha))})
    tmp = candidate[order(pel)[1]]
    inform = c(inform, tmp)
  }
  return(inform)
}


# #Test
# dat = as.matrix(read.table("subdata1.txt"))
# a = rep(mean(dat),dim(dat)[2])
# b = 1-a
# warm = 500
# iter = 100
# step = 10
# size = 100
# alpha = 15
# system.time({
#   cl_sample = dpmm_beta(a,b,x0 = dat, warm=warm, size=size,iter=iter, 
#                         step=step, alpha = alpha)
# })
# system.time({
#   inform = inform_beta(cl_sample,iter,size, nA = 4, nT = 5,a,b, x0 = dat,alpha)
# })
