#source("dpmm_beta.R")
#source("Update_beta.R")
#source("Initial_beta.R")
source("pel1_beta.R")
source("pel2.R")
source("inform_beta_v1.R")
source("inform_beta_v2.R")
#load data
load("cl_GDSC.RData")



#choose iterations, warm up step, step length
warm = 500
iter = 100
step = 10
size = 2000

#find an empirical best prior mass alpha
alpha = 3
#alpha = 5
# Leave-one-out cross validation for BOISE framework
### Lets try nA = 5, nT = 20
nA = 7
nT = 20

value <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(value) + 1

test = GDSC[i, ]
train = GDSC[-i, ]
a = rep(mean(train),dim(train)[2])
b = 1 - a

### Sample for xi

#cl_sample = dpmm_beta(a,b,x0 = train, warm, size=1, iter, step, alpha)
cl_sample = cl_GDSC[[i]]
rm(cl_GDSC)
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


tmp = read.table("GDSC5_result.txt")
pre_inform = as.numeric(unlist(strsplit(as.character(tmp$V2[i]), split = ' ')))

#inform = inform_beta1(cl_sample, P, iter,size, nA = nA, nT = nT,a,b, x0 = train,alpha)
inform = inform_beta2(cl_sample,P, iter,size,nT,a,b,x0 = train,alpha,
                      inform = pre_inform, nAdd = 2)

inform = paste(inform, collapse = " ")
result = data.frame("A" = i, "Inform" = inform) 

## Write results
if(nchar(value) == 1){
  write.table(result,file = paste("GDSC7_result00", value, ".txt", sep = ""),col.names = F,row.names = F)
} else if(nchar(value) == 2){
  write.table(result,file = paste("GDSC7_result0", value, ".txt", sep = ""),col.names = F,row.names = F)
} else{
  write.table(result,file = paste("GDSC7_result", value, ".txt", sep = ""),col.names = F,row.names = F)
}
# svname <- paste("clust_tar",i,sep="_")
# svname <- paste(svname, "RData", sep=".")
# save( i, cl_sample, file=svname)
