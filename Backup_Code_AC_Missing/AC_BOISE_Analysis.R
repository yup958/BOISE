source("Initial_beta.R")
source("Update_beta.R")
source("dpmm_beta.R")
source("clust_sum.R")
source("pel1_beta.R")
source("pel2.R")
source("inform_beta_v1.R")
source("inform_beta_v2.R")
source("Entropy.R")
source("Evaluate.R")
#load data
load("pkis1.rda")
## Use 2sd criteria to create binary matrix
dat <- t(apply(pkis1, 1, function(x){
  thres = mean(x) + 2 * sd(x)
  return(as.numeric(x>thres))
}))
rm(pkis1)

#choose iterations, warm up step, step length
warm = 500
iter = 100
step = 10
size = 1000
#find an empirical best prior mass alpha
alpha = 15

# One experiment in leave-one-out cross validation for BOISE framework
# First create cl_sample objects for ith training set
value <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(value) + 1
test = dat[i, ]
train = dat[-i, ]
a = rep(mean(train),ncol(train))
b = 1 - a
m = ncol(train)
n = nrow(train)
cl_sample = dpmm_beta(a,b,x0 = train, warm, iter, step, alpha)

### Create a matrix summarize information of each clustering.
### K by (n+1) matrix. Each row for one cluster, last column as number of targets in cluster
P = clust_sum(cl_sample, train, iter, a, b)

### Sample possible x_i*
XX = rep(0, iter*size*m)
dim(XX) = c(iter, size, m)
for (j in 1:iter) {
  K = cl_sample$KK[j]
  p = rep(0, K + 1)
  p[K + 1] = alpha / (n + alpha)
  p[1:K] = P[[j]][,m+1] / (n+alpha)
  XX[j, , ] = t(as.matrix(sapply(1:size, function(s){
    classi = which(rmultinom(1, 1, p) == 1)
    if(classi == K+1){
      post_theta = a / (a + b)
    } else{
      post_theta = P[[j]][classi,1:m]
    }
    new_xi = sapply(1:m, function(x){
      return(as.numeric(rbinom(1, 1, p = post_theta[x])))
    })
    return(new_xi)
  })))
}

## Accelerated BOISE
candidate = 1:ncol(train)
info_gain = rep(0,length(candidate))
info_gain= sapply(candidate, function(j){
  inform = j
  tmp = sapply(1:iter, function(k){
    Pk = P[[k]]
    XXk = XX[k, , ]
    return(Entropy(P = Pk,inform = inform,XX = XXk,alpha,a,b))
  })
  return(mean(tmp))
})
inform = order(info_gain)[1]
# tmp = read.table("PKIS18_result.txt")
# inform = as.numeric(unlist(strsplit(as.character(tmp$V2[i]), split = ' ')))
candidate = candidate[-inform]

while((length(inform) < nA)){
  info_gain = rep(0,length(candidate))
  info_gain= sapply(candidate, function(j){
    inform = c(inform,j)
    tmp = sapply(1:iter, function(k){
      Pk = P[[k]]
      XXk = XX[k, , ]
      return(Entropy(P = Pk,inform = inform,XX = XXk,alpha,a,b))
    })
    return(mean(tmp))
  })
  best = candidate[order(info_gain)[1]]
  candidate = candidate[-order(info_gain)[1]]
  inform = c(inform,best)
}

inform = paste(inform, collapse = " ")
result = data.frame("A" = i, "Inform" = inform) 

## Write results
if(nchar(value) == 1){
  write.table(result,file = paste("PKIS1_16_result00", value, ".txt", sep = ""),col.names = F,row.names = F)
} else if(nchar(value) == 2){
  write.table(result,file = paste("PKIS1_16_result0", value, ".txt", sep = ""),col.names = F,row.names = F)
} else{
  write.table(result,file = paste("PKIS1_16_result", value, ".txt", sep = ""),col.names = F,row.names = F)
}