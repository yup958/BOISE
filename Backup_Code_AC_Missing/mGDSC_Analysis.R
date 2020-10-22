rm(list = ls())
## data cleaning
GDSC = read.csv("GDSC.csv")
rownames(GDSC) = GDSC$X
GDSC$X = NULL
GDSC = apply(GDSC, 2, function(x){
  return(1*(x<=-2))
})

missing = apply(GDSC,1, function(x){
  return(sum(is.na(x)))
})
test.index = which(missing == 0)
test = GDSC[test.index,]
train = GDSC[-test.index,]

# beta prior
a = rep(mean(train,na.rm = T),ncol(train))
b = 1-a
#choose iterations, warm up step, step length
warm = 500
iter = 100
step = 10
size = 1000

#find an empirical best prior mass alpha
alpha = 15

# Leave-one-out cross validation for BOISE framework
### Lets try nA = 1,2,3, nT = 10
nA = 8
nT = 31

# test dpmm
source("Initial_beta.R")
source("Update_beta.R")
source("dpmm_beta.R")
cl_sample = dpmm_beta(a,b,train,warm,iter, step, alpha)

source("clust_sum.R")
P = clust_sum(cl_sample, train, iter, a, b)
## Sample for x_i*
m = ncol(train)
n = nrow(train)
cl_sample$XX = rep(0, iter*size*m)
dim(cl_sample$XX) = c(iter, size, m)
for (j in 1:iter) {
  K = cl_sample$KK[j]
  p = rep(0, K + 1)
  p[K + 1] = alpha / (n + alpha)
  p[1:K] = P[[j]][,m+1] / (n+alpha)
  cl_sample$XX[j, , ] = t(as.matrix(sapply(1:size, function(s){
    classi = which(rmultinom(1, 1, p) == 1)
    if(classi == K+1){
      post_theta = a / (a + b)
    } else{
      post_theta = P[[j]][classi,1:m]
    }
    new_xi = rbinom(m,1,post_theta)
    return(new_xi)
  })))
}

source("pel1_beta.R")
source("pel2.R")
iter = 100
size = 1000
nA = 8
nT = 31


#New cluster method no imputation:
#inform = c(100, 299,  96, 247,  93, 110,  40,  22)
inform = inform_beta1(cl_sample,P,iter,size,nA,nT,a,b,x0=train,alpha)
inform = 100
candidate = (1:ncol(train))[-inform]
n = nrow(train)
m = ncol(train)
cl_sample$XX = rep(0, iter*size*m)
cl_sample$XX[actives] = 1
dim(cl_sample$XX) = c(iter, size, m)

# Distributed computing for BOISE:
pel1 = unlist(mclapply(candidate, function(k){
  tmp = pel1_beta(cl_sample, P, iter, size, A = c(inform,k), 
                  nA = 1+length(inform), nT, a, b, x0 = train, alpha)
  return(tmp)
},mc.cores = detectCores()))
inform = c(inform, candidate[order(pel1)[1]])

#Evaluation
source("Evaluate.R")
nef.result = rep(0,23)
nef.result = unlist(mclapply(1:23,function(i){
  return(Evaluate(P,inform,"nef",test[i,],train,
                  nA=length(inform),nT,iter,a,b,alpha))
},mc.cores = detectCores()))

# Distributed computing for AC-BOISE:
source("~/RAwork/esdd/BOISE/Code/Entropy2.R")
candidate = 1:ncol(train_imp)
info_gain = rep(0,length(candidate))
info_gain = unlist(mclapply(candidate, function(j){
  inform = j
  tmp = 0
  for (k in 1:iter) {
    tmp_cl = list(K=cl_sample$KK[k],N = cl_sample$NN[k,1:cl_sample$KK[k]],
                  C = cl_sample$CC[k,])
    tmp = tmp + Entropy2(cl=tmp_cl,inform = inform,dat = train_imp)
  }
  return(tmp/iter)
}, mc.cores = detectCores()))

inform = order(info_gain)[1]
candidate = candidate[-inform]

while((length(inform) < nA)){
  info_gain = rep(0,length(candidate))
  info_gain = unlist(mclapply(candidate, function(j){
    inform = c(inform,j)
    tmp = 0
    for (k in 1:iter) {
      tmp_cl = list(K=cl_sample$KK[k],N = cl_sample$NN[k,1:cl_sample$KK[k]],
                    C = cl_sample$CC[k,])
      tmp = tmp + Entropy2(cl=tmp_cl,inform = inform,dat = train_imp)
    }
    return(tmp/iter)
  }, mc.cores = detectCores()))
  best = candidate[order(info_gain)[1]]
  candidate = candidate[-order(info_gain)[1]]
  inform = c(inform,best)
}
## AC-BOISE inform = c(100,2,8,22,284,247,1,195)
# BF_w
nef.result = rep(0,23)
auc.result = nef.result
nA=8
nT=31
dat = read.csv("~/JupyterNote/GDSC_knn.csv", header = F)
GDSC = read.csv("~/Drug_Discovery/BOISE/CR_BetaPrior/GDSC.csv")
rownames(GDSC) = GDSC$X
GDSC$X = NULL
rownames(dat) = rownames(GDSC)
colnames(dat) = colnames(GDSC)
missing = apply(GDSC,1, function(x){
  return(sum(is.na(x)))
})
test.index = which(missing == 0)

train = dat[-test.index,]
test = dat[test.index,]
test = apply(test,2,function(x){
  return(1*(x<=-2))
})
rm(GDSC)
hits = apply(train, 2, sum)
inform = order(hits,decreasing = T)[1:nA]
for (i in 1:23) {
  xA = test[i, inform]
  Score = apply(train, 2,function(x){
    corr = cor(x, train[,inform])
    return(sum(corr * xA))
  })
  Score[inform[which(xA==1)]] = rep(max(Score) + 1, sum(xA))
  Score[inform[which(xA==0)]] = rep(min(Score) - 1, nA - sum(xA))
  Score = as.vector(Score)
  rocobj = roc(test[i,],Score)
  auc.result[i] = auc(rocobj)
  tophit = sum(test[i,order(Score,decreasing = T)[1:nT]])
  hit = sum(test[i,])
  maxhit = min(hit,nT)
  nef.result[i] = ((tophit/nT - hit/ncol(train)) / (maxhit/nT - hit/ncol(train)) + 1)/2
  
}

## RS result
x =read.table("~/Documents/MATLAB/Regression Selection/2_NEF_result.txt",sep = ",")
x = as.matrix(x)
nef.result = rep(0,23)
auc.result = nef.result
for (i in 1:23) {
  auc.result[i] = pROC::roc( test[i,],x[i,])$auc
  top = order(x[i,])[1:31]
  predict.hit = sum(test[i,top])
  hit = sum(test[i,])
  max.hit = max(hit, 31)
  nef.result[i] = ((predict.hit/31 - hit/304)/(max.hit/31 - hit/304)+1)/2
}

