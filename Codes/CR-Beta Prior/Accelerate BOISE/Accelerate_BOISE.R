setwd("~/Drug_Discovery/BOISE/CR_BetaPrior/")
rm(list=ls())
load("cl.RData")
load("newtar.RData")
load("cl_GDSC.RData")
source("pel1_beta.R")
source("pel2.R")
library(pROC)
dat = GDSC
auc.result = rep(0,nrow(dat))
nef.result = rep(0,nrow(dat))
nA = 8
nT = 20

for (i in 1:nrow(dat)) {
  print(i)
  info_gain = rep(0,ncol(dat))
  train = dat[-i,]
  test = dat[i,]
  cl_sample = cl_GDSC[[i]]
  n = nrow(train)
  
  for (k in 1:iter) {
    tmp_cl = list(K=cl_sample$KK[k],N = cl_sample$NN[k,1:cl_sample$KK[k]], C = cl_sample$CC[k,])
    p = tmp_cl$N / n
    Entropy0 = sum(-p*log(p))
    Entropy1 = sapply(1:ncol(train), function(j){
      aj = sapply(1:tmp_cl$K, function(x){
        target = which(tmp_cl$C == x)
        return(sum(train[target,j]))
      })
      bj = tmp_cl$N - aj
      pos_p = aj / sum(aj)
      pos_p = pos_p[which(pos_p > 0)]
      pos_entropy = (sum(aj)/n) * sum(-pos_p*log(pos_p))
      neg_p = bj / sum(bj)
      neg_p = neg_p[which(neg_p>0)]
      neg_entropy = (sum(bj)/n) * sum(-neg_p*log(neg_p))
      return(pos_entropy+neg_entropy)
    })
    info_gain = info_gain + Entropy0 - Entropy1
  }
  
  candidate = order(info_gain,decreasing = T)
  inform = candidate[1]
  point = 2
  while((length(inform) < nA) && (point < length(candidate))){
    correlation = cor(train[,inform],train[,candidate[point]])
    inactive = 1*(correlation > 0.2)
    if(sum(inactive)>0){
      point = point + 1
    } else{
      inform = c(inform, candidate[point])
      point = point + 1
    }
  }
  
  a = rep(mean(train),dim(train)[2])
  b = 1 - a
  Score = rep(0, dim(train)[2])
  xA = test[inform]
  for (k in 1:iter) {
    tmp_cl = list(K = cl_sample$KK[k], N = cl_sample$NN[k,], C = cl_sample$CC[k,])
    post_theta = pel2_beta(tmp_cl, x0 = train, xA = test[inform], nA = 16, A = inform, nT = 36, a, b, alpha=15)
    Score = Score + post_theta
  }
  Score[inform[which(xA==1)]] = rep(max(Score) + 1, sum(xA))
  Score[inform[which(xA==0)]] = rep(min(Score) - 1, nA - sum(xA))
  top = order(Score,decreasing = T)[1:nT]
  result = sum(test[top])
  
  hit = sum(test)
  maxhit = min(hit,nT)
  nef.result[i] = ((result/nT - hit/207) / (maxhit/nT - hit/207) + 1)/2
  
  Score = as.vector(Score)
  rocobj = roc(test,Score)
  auc.result[i] = rocobj$auc
}
testresult_nef = read.csv("pkis18_eval_NEF10.csv")
testresult_roc = read.csv("pkis18_eval_ROCAUC.csv")

testresult_nef = read.csv("~/Drug_Discovery/BOISE/CR_BetaPrior/CR_BetaPrior_Result/pkis1_eval_NEF10.csv")
testresult_roc = read.csv("~/Drug_Discovery/BOISE/CR_BetaPrior/CR_BetaPrior_Result/pkis1_eval_ROCAUC.csv")

testresult_nef = read.csv("~/Drug_Discovery/BOISE/CR_BetaPrior/GDSC8_eval_NEF10.csv")
testresult_roc = read.csv("~/Drug_Discovery/BOISE/CR_BetaPrior/GDSC8_eval_ROCAUC.csv")

write.csv(nef.result,file = "nef8.csv",row.names = F)
write.csv(auc.result,file = "roc8.csv",row.names = F)
