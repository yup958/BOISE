### R script submitted to CHTC for HPC computing
### "newtar.RData" is a sampled CR clustering assignments with size 100, on whole PKIS1 data.
### "inform" is informer set obtained in last step (in this script, an informer set of size 13)
### An initial informer set could be selected by inform_beta1 function
### This script will return the pel1 loss on all candidate compounds, 
### then we will choose the next informer compound with least pel1 loss.
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
nA = 14
nT = 36

test = dat[225:227,]
train = dat[1:224,]
inform = c(217,252,8,33,42,178,67,88,35,162,34,126,180,278)
a = rep(mean(train),dim(train)[2])
b = 1 - a
#inform = inform_beta1(cl_sample,iter,nA = nA, nT = nT,a,b, x0 = train,alpha)
#inform = inform_beta2(cl_sample,iter,nT,a,b,x0 = train,alpha,
                      inform = inform, nAdd = 1)
candidate = (1:dim(train)[2])[-inform]
value <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(value) + 1
candidate = candidate[i]
step = 1
n = length(inform)
pel = pel1_beta(cl_sample,iter,A=c(inform,candidate),nA=n+step,nT,a,b,x0=train,alpha)
result = data.frame("A" = candidate, "pel1" = pel) 

if(nchar(value) == 1){
  write.table(result,file = paste("pel_result00", value, ".txt", sep = ""),col.names = F,row.names = F)
} else if(nchar(value) == 2){
  write.table(result,file = paste("pel_result0", value, ".txt", sep = ""),col.names = F,row.names = F)
} else{
  write.table(result,file = paste("pel_result", value, ".txt", sep = ""),col.names = F,row.names = F)
}
