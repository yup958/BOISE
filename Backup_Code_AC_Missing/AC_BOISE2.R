### Accelerate BOISE selection for initial test on data sets
## Input: sampled DPMM cluster assginments cl_sample;
##        informer set size nA;
##        training set train
## Output: informer set inform of size nA.

AC_BOISE2 <- function(cl_sample, nA, train){
  source("Entropy2.R")
  candidate = 1:ncol(train)
  info_gain = rep(0,length(candidate))
  info_gain= sapply(candidate, function(j){
    inform = j
    tmp = sapply(1:iter, function(k){
      tmp_cl = list(K=cl_sample$KK[k],N = cl_sample$NN[k,1:cl_sample$KK[k]], C = cl_sample$CC[k,])
      return(Entropy2(cl=tmp_cl,inform = inform,dat = train))
    })
    return(mean(tmp))
  })
  inform = order(info_gain)[1]
  candidate = candidate[-inform]
  while((length(inform) < nA)){
    info_gain = rep(0,length(candidate))
    info_gain= sapply(candidate, function(j){
      inform = c(inform,j)
      tmp = sapply(1:iter, function(k){
        tmp_cl = list(K=cl_sample$KK[k],N = cl_sample$NN[k,1:cl_sample$KK[k]], C = cl_sample$CC[k,])
        return(Entropy2(cl=tmp_cl,inform = inform,dat = train))
      })
      return(mean(tmp))
    })
    best = candidate[order(info_gain)[1]]
    candidate = candidate[-order(info_gain)[1]]
    inform = c(inform,best)
  }
  return(inform)
}

# nA = 16
# system.time({
#   inform =  AC_BOISE2(cl_sample, nA, train)
#   })
# 12537
