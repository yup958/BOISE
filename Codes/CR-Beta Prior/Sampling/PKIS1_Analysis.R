### R script submitted to CHTC
### cl.RData is a list of 100 sampled clustering assignments on each of 224 leave-one-out cross validation experiment. 
### That is, a list of length 224, each element is also a list with 3 elements {KK,NN,CC}
### Will return a larger informer set for each of the leave-one-out experiment.
source("pel1_beta.R")
source("pel2.R")
source("inform_beta_v1.R")
source("inform_beta_v2.R")
#load data
# cl = list(rep(0,224))
# for (i in 1:224) {
#   test = dat[i, ]
#   train = dat[-i, ]
#   a = rep(mean(train),dim(train)[2])
#   b = 1 - a
#   cl_sample = dpmm_beta(a,b,x0 = train, warm=500, size=1, iter=100, step=10, alpha=15)
#   cl[[i]] = cl_sample
# }
load("cl.RData")
## Use 2sd criteria to create binary matrix
# dat <- t(apply(pkis1, 1, function(x){
#   thres = mean(x) + 2 * sd(x)
#   return(as.numeric(x>thres))
# }))
# dat = dat[ ,-which(apply(dat,2,sum) == 0)]
# rm(pkis1)

#choose iterations, warm up step, step length
warm = 500
iter = 100
step = 10
size = 1000

#find an empirical best prior mass alpha
alpha = 15
# Leave-one-out cross validation for BOISE framework
nA = 7
nT = 36

value <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(value) + 1

test = dat[i, ]
train = dat[-i, ]
a = rep(mean(train),dim(train)[2])
b = 1 - a

### Sample for xi, possible outcomes on all 366 compounds

#cl_sample = dpmm_beta(a,b,x0 = train, warm, size, iter, step, alpha)
cl_sample = cl[[i]]
n = dim(train)[1]
m = dim(train)[2]
cl_sample$XX = rep(0, iter*size*m)
dim(cl_sample$XX) = c(iter, size, m)
for (k in 1:iter) {
  p = rep(0, cl_sample$K[k] + 1)
  p[cl_sample$K[k] + 1] = alpha / (n + alpha)
  p[1:cl_sample$K[k]] = sapply(1:cl_sample$K[k],function(x){
    return(cl_sample$N[k,x]/(n+alpha))})
  cl_sample$XX[k, , ] = t(as.matrix(sapply(1:size, function(s){
    classi = which(rmultinom(1, 1, p) == 1)
    if(classi == cl_sample$K[k]+1){
      posta = a
      postb = b
    } else{
      targets = which(cl_sample$C[k,] == classi)
      posta = a
      postb = b
      for (t in targets) {
        posta = posta + train[t, ]
        postb = postb + 1 - train[t,]
      }
    }
    p = sapply(1:m, function(x){
      return(rbeta(1,posta[x],postb[x]))
    })
    new_xi = sapply(1:m, function(x){
      return(as.numeric(rbinom(1, 1, p = p[x])))
    })
    return(new_xi)
  })))
}


tmp = read.table("Test6_result.txt")
pre_inform = as.numeric(unlist(strsplit(as.character(tmp$V2[i]), split = ' ')))

#inform = inform_beta1(cl_sample,iter,size, nA = nA, nT = nT,a,b, x0 = train,alpha)
inform = inform_beta2(cl_sample,iter,size,nT,a,b,x0 = train,alpha,
                      inform = pre_inform, nAdd = 1)
# Score = rep(0, dim(train)[2])
# for (k in 1:iter) {
#   tmp_cl = list(K = cl_sample$KK[k], N = cl_sample$NN[k,], C = cl_sample$CC[k,])
#   post_theta = pel2_beta(tmp_cl, x0 = train, xA = test[inform], nA = nA, A = inform, nT = nT, a, b, alpha)
#   Score = Score + post_theta
# }
# top = order(Score,decreasing = T)[1:nT]
# result = data.frame("A" = i, "Accuracy" = sum(test[top]))
inform = paste(inform, collapse = " ")
result = data.frame("A" = i, "Inform" = inform) 

## Write results
if(nchar(value) == 1){
  write.table(result,file = paste("Test7_result00", value, ".txt", sep = ""),col.names = F,row.names = F)
} else if(nchar(value) == 2){
  write.table(result,file = paste("Test7_result0", value, ".txt", sep = ""),col.names = F,row.names = F)
} else{
  write.table(result,file = paste("Test7_result", value, ".txt", sep = ""),col.names = F,row.names = F)
}
