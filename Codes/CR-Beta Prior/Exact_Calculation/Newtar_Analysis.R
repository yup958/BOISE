library(R.utils)
source("pel1_beta.R")
source("pel2.R")
source("inform_beta_v1.R")
source("inform_beta_v2.R")
load("newtar.RData")

warm = 500
iter = 100
step = 10
size = 1000
alpha = 15
nA = 7
nT = 36

test = dat[225:227,]
train = dat[1:224,]
a = rep(mean(train),dim(train)[2])
b = 1 - a
#inform = inform_beta1(cl_sample,iter,nA = nA, nT = nT,a,b, x0 = train,alpha)
inform = inform_beta2(cl_sample,iter,nT,a,b,x0 = train,alpha,
                      inform = inform, nAdd = 1)
result = data.frame("A" = 1, "Inform" = inform) 
write.table(result,file = "Test7_result.txt",col.names = F,row.names = F)