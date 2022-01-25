setwd('~/RAwork/BOISE/BOISE_followup/FDA_analysis/')
load('fda_data_rearranged.RData')
set.seed(2)
library(BOISE)
block_res$row_block = block_res$row_block[1:5,]
block_res$col_block = block_res$col_block[1:5,]
active_ind = apply(block_res$col_block, 2, sum)
active_ind = which(active_ind > 0)
cpds = colnames(block_res$col_block[, active_ind])
dat = dat[, cpds]

## Train-test separation
count = apply(dat,1, function(x){return(sum(is.na(x)))})
test_AID = names(sort(count)[1:100])
train_AID = setdiff(rownames(dat), test_AID)
train = dat[train_AID, ]
test = dat[test_AID,]
rm(dat)

m0 = 10
a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
sample_size = 500
system.time({
  cl_sample = dpmm_beta(train, a, b, m0= m0, burn_in = 1000, sample_size=sample_size, thinning = 15)
})

## Baseline ROCAUC, no informers, just score with hit frequency
baseline_roc = rep(NA,100)
for (i in 1:100){
  Score = apply(train,2,function(x){return(mean(x, na.rm = T))})
  complete_idx = which(!is.na(test[i,]))
  tmp = as.vector(test[i,complete_idx])
  Score = as.vector(Score[complete_idx])
  if(sum(test[i,], na.rm = T) == 0){
    print('All observed responses are 0!')
  } else{
    rocobj = pROC::roc(response = tmp, predictor = Score, quiet = TRUE)
    baseline_roc[i] = rocobj$auc
  }
}
summary(baseline_roc)

### load clustering results, informer selection from Condor
load('clust_res_10.RData')
inform = c(25, 195, 264, 285, 525, 189, 102, 259, 32, 205, 112, 71,103, 95, 138,141)
inform_CID = colnames(train)[inform]
nT = as.integer(ncol(train) * 0.1)

### Test set performance
## revise Evaluate function to incorporate NAs
Evaluate <-
  function(cl_sample, inform, measure, percent,
           test, train, nT, sample_size, alpha, beta, m0){
    P = clust_sum(cl_sample,train,sample_size, alpha, beta)
    Score = rep(0, ncol(train))
    complete_idx = which(!is.na(test))
    inform = intersect(inform, complete_idx)
    xA = test[inform]
    nA=length(inform)
    m = ncol(train)
    post_probs = matrix(0, 1, sample_size)
    post_thetas = matrix(0, sample_size, m)
    for (k in 1:sample_size){
      postls = pel2_beta(P[[k]], x0=train, xA, A=inform, nT, alpha, beta, m0)
      post_probs[k] = postls$post_prob
      post_thetas[k, ] = postls$post_theta
    }
    post_probs = post_probs / (sum(post_probs))
    Score = post_probs %*% post_thetas
    Score[inform[which(xA==1)]] = rep(max(Score) + 1, sum(xA, na.rm = T))
    Score[inform[which(xA==0)]] = rep(min(Score) - 1, sum(1-xA, na.rm = T))
    test = as.vector(test[complete_idx])
    Score = as.vector(Score[complete_idx])
    if(measure == "nef"){
      if(sum(test, na.rm = T) == 0){
        print('All observed responses are 0!')
        result = NA
      } else{
        nTop = round(ncol(train) * percent, 0)
        top = order(Score,decreasing = T)[1:nTop]
        pred_hit = sum(test[top])
        hit = sum(test)
        maxhit = min(hit,nTop)
        result = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
      }
    } else if(measure == 'selectivity_nef'){
      nTop = round(ncol(train) * percent, 0)
      top = order(Score,decreasing = T)[1:nTop]
      selectivity_score = colMeans(train, na.rm = T)
      test = test - selectivity_score
      pred_hit = sum(test[top])
      hit = sum(test)
      maxhit = sum(sort(test, decreasing = T)[1:nTop])
      result = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
    } else if(measure == "rocauc"){
      if(sum(test, na.rm = T) == 0){
        print('All observed responses are 0!')
        result = NA
      } else{
        rocobj = pROC::roc(response = test, predictor = Score, quiet = TRUE)
        result = rocobj$auc
      }
    } else if(measure %in% c("mat", "f")){
      pred.obj = ROCR::prediction(Score, test)
      perform.obj = ROCR::performance(pred.obj, measure)
      result = max(unlist(perform.obj@y.values),na.rm = T)
    } else{
      print("Criteria is not supported.")
      result = 0
    }
    return(result)
  }
roc_results = read.table('roc_results.txt', sep = ' ', header=T)
nA_results = read.table('nA_results.txt', sep = ' ', header=T)
nef_results = read.table('nef_results.txt', sep = ' ', header=T)

for (k in 16:16) {
  nA_name = paste('nA_', as.character(k), sep = '')
  roc_name = paste('rocauc_', as.character(k), sep = '')
  #nef_name = paste('nef_', as.character(k), sep = '')
  roc_results[, roc_name] = rep(NA, 100)
  #nef_results[, nef_name] = rep(NA, 100)
  nA_results[, nA_name] = rep(0, 100)
  tmp_inform = inform[1:k]
  for (i in 1:100) {
    complete_idx = which(!is.na(test[i,]))
    valid_inform = intersect(complete_idx, tmp_inform)
    nA_results[i, nA_name] = length(valid_inform)
    if (length(valid_inform) >= 1){
      roc_results[i, roc_name] = Evaluate(cl_sample, tmp_inform, measure = 'rocauc', percent = 0.1,
                                       test=test[i,], train, nT, sample_size, a, b, m0)
      #nef_results[i, nef_name] = Evaluate(cl_sample, tmp_inform, measure = 'nef', percent = 0.1,
      #                                    test=test[i,], train, nT, sample_size, a, b, m0)
    }
  }
}
write.table(nA_results, file = 'nA_results.txt',row.names = F)
write.table(roc_results, file = 'roc_results.txt',row.names = F)
write.table(nef_results, file = 'nef_results.txt',row.names = F)

## drop from 6 to 7
worse_set = which(roc_results$rocauc_6 - roc_results$rocauc_7 > 0.05)
worse_AID = rownames(test)[worse_set]
test[worse_AID, inform_CID[1:7]]
block_res$col_block[, inform_CID[1:7]]
### take some clustering for example
grp = 10
for (k in 1:cl_sample$KK[grp]) {
  grp_id = which(cl_sample$CC[grp,] == k)
  size = cl_sample$NN[grp,k]
  avg_rate = mean(train[grp_id, inform_CID[7]], na.rm = T)
  print(paste('Clust:', k, 'Size:', size, 'Avg_rate:', avg_rate, sep=' '))
}
### average positive in test/train set
summary(test[,inform_CID[7]])
summary(train[,inform_CID[7]])

### Plots
means = apply(roc_results, 2,function(x){return(mean(x, na.rm = T))})
medians = apply(roc_results, 2,function(x){return(median(x, na.rm = T))})
par(mfrow = c(2,1))
plot(y = means[3:15], x=2:14, xlab = 'Informer_size', ylab = 'ROCAUC mean')
lines(y=means[3:15], x=2:14)
abline(h = mean(baseline_roc, na.rm = T), col = 'red')
plot(y=medians[3:15], x=2:14, xlab = 'Informer_size', ylab = 'ROCAUC median')
lines(y=medians[3:15], x=2:14)
abline(h = median(baseline_roc, na.rm = T), col = 'red')

### further investigation: random sampled informer set, will same thing happen?
roc_results = read.table('roc_results.txt', sep = ' ', header=T)
nA_results = read.table('nA_results.txt', sep = ' ', header=T)
nsamp = 1
#prob_weight = apply(train,2, function(x){return(sum(!is.na(x)))})
prob_weight = apply(train,2, function(x){return(mean(x, na.rm = T))})

for (s in 1:nsamp) {
  inform = sample(1:933, 14, prob = prob_weight)
  for (k in 3:14) {
    tmp_inform = inform[1:k]
    nA_name = paste('nA_', as.character(k), sep = '')
    roc_name = paste('rocauc_', as.character(k), sep = '')
    roc_results[, roc_name] = rep(NA, 100)
    nA_results[, nA_name] = rep(0, 100)
    for (i in 1:100) {
      complete_idx = which(!is.na(test[i,]))
      valid_inform = intersect(complete_idx, tmp_inform)
      nA_results[i, nA_name] = nA_results[i, nA_name] + length(valid_inform)
      if (length(valid_inform) >= 1){
          roc_results[i, roc_name] = Evaluate(cl_sample, tmp_inform, measure = 'rocauc', percent = 0.1,
                                              test=test[i,], train, nT, sample_size, a, b, m0) 
          }
      }
    }
}

