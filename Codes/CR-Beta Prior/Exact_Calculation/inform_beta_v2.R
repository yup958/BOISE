### Informer set selection function version 2
### Pre-step: Sample DPMM samples of cl_sample. 
###           A previous selected smaller informer set
### Input: list cl_sample, corresponding parameters a,b,iter,size,alpha; nT; x0
###         Pre-selected informer set inform, number of informers to add nAdd
### Output: Informer set A with size n+nAdd (via advanced adaptive selection)

inform_beta2 <- function(cl_sample, iter, nT, a, b, x0, alpha, inform, nAdd){
  source("pel1_beta.R")
  step = 1
  n = length(inform)
  candidate = (1:dim(x0)[2])[-inform]
  while (step <= nAdd) {
    pel = rep(0,length(candidate))
    pel = sapply(candidate, function(x){
      return(pel1_beta(cl_sample, iter, A = c(inform,x), nA = n+step,nT,a,b,x0, alpha))})
    tmp = candidate[order(pel)[1]]
    inform = c(inform, tmp)
    step = step + 1
    candidate = candidate[-which(candidate == tmp)]
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
