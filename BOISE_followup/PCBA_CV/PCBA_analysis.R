set.seed(817)
setwd('~/RAwork/BOISE/BOISE_followup/PCBA_CV/')
dat = as.matrix(data.table::fread('~/Downloads/condensed_subset_PCBA128.csv.gz'),rownames=1)
dat = dat[, !colnames(dat) %in% c("rdkit SMILES","Morgan FP_2_1024","Cluster ID")]
col_name = colnames(dat)
row_name = rownames(dat)
dat = matrix(as.numeric(dat), ncol = ncol(dat))
colnames(dat) = col_name
rownames(dat) = row_name
dat = t(dat)
### overlaps with FDA
load('~/RAwork/BOISE/BOISE_followup/FDA_CV/fda_data.RData')
fda_cids = colnames(dat)
fda_cids = sapply(fda_cids, function(s){return(unlist(strsplit(s, '_'))[2])})
common_cids = intersect(fda_cids, row_name)
rm(dat)

counts = colSums(dat, na.rm = T)
dat = dat[, which(counts > 0)]
cids = colnames(dat)

### chemical clustering
dat = as.matrix(data.table::fread('~/Downloads/condensed_subset_PCBA128.csv.gz'),rownames=1)
dat = dat[rownames(dat)%in%cids,]
fps = dat[,'Morgan FP_2_1024']
rm(dat)
fps = t(sapply(fps, function(x){return(as.numeric(unlist(strsplit(x,''))))}))
D = dist(fps, method = 'binary')
load('PCBA_dist.RData')
cl = cluster::diana(x = D, diss = T, stop.at.k = 20, keep.diss = F)
cl = hclust(D, method = 'average')
silhouette_score <- function(k){
  tmp_cl = cutree(cl, k)
  ss <- cluster::silhouette(x = tmp_cl, dist = cl$diss)
  mean(ss[, 3])
}
k <- 2:25
avg_sil <- sapply(k, silhouette_score)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
tmp_cl = cutree(cl, k=23)
table(tmp_cl)
res = data.frame(cid = colnames(dat), clust = tmp_cl)
write.csv(res, file = 'chemical_clustering_res.csv', row.names = F)

### nA = 1000
### idx = 1: NEF10 0.58;
### idx = 2: NEF10
# function to compute validity score
valid_score <- function(k, D=D){
  kmed = cluster::pam(D, k, diss = T)
  # among cluster distance
  inter = 0
  for(i in 1:(k-1)){
    medoid1 = kmed$medoids[i]
    for (j in (i+1):k) {
      medoid2 = kmed$medoids[j]
      inter = inter + D[medoid1, medoid2]
    }
  }
  inter = inter / (k*(k-1) / 2)
  # within cluster distance
  m = nrow(D)
  intra = rep(0, m)
  for (i in 1:m) {
    label = kmed$clustering[i]
    med = kmed$medoids[label]
    intra[i] = D[i, med]
  }
  intra = mean(intra)
  return(intra /inter)
}
harmonic_sum <- function(start, end){
  res = 0
  for(k in start:end){
    res = res + 1/k
  }
  return(res)
}

idxs = sample(1:102)#102
sample_size = 500
for(idx in idxs[31:31]){
  print(idx)
  test_id = rownames(dat)[idx]
  test = dat[rownames(dat)==test_id,]
  complete_idx = which(!is.na(test))
  train = dat[!row.names(dat)==test_id, complete_idx]
  test = test[complete_idx]
  ### Approximate clustering with random distance using Hamming / Manhatton dist or Jaccard dist
  m = nrow(train)
  D = matrix(0,m,m)
  for (i in 1:m) {
    for (j in 1:m) {
      #D[i, j] = sum(abs(train[i,] - train[j,]), na.rm = T) # Hamming/Manhatton
      set1 = which(train[i,] == 1)
      set2 = which(train[j,] == 1)
      if(length(union(set1, set2)) == 0)
        D[i, j] = 0
      else
        D[i, j] = 1 - length(intersect(set1, set2)) / length(union(set1, set2))
    }
  }
  ## randomized distances
  # empirically estimate a0 and a1
  dists = as.vector(D)
  dists = dists[-which(dists == 0)]
  d0 = mean(1/dists) / var(1/dists)
  log_obj <- function(a){
    m = 101
    const = -lbeta(a[1], a[2])
    res = 0
    for (i in 1:(m-1)) {
      for (j in (i+1):m){
        if(D[i, j] > 0)
          res = res + const + a[1]*log(d0) + (a[2]-1)*log(D[i, j]) + 
            a[2]*log(a[2]) - (a[1]+a[2])*log(d0+a[2]*D[i,j])
      }
    }
    return(-res)
  }
  optim_obj = nlminb(start = c(250,250), objective = log_obj, lower = 0, upper = 1e6,control = list(eval.max = 600))
  a = sum(optim_obj$par) #3.27e5 for Hamming, 2.95e5 for Jaccard
  
  CC = matrix(0, sample_size, nrow(train))
  for (id in 1:sample_size) {
    epsilons = rgamma(n=m, shape = a/2, rate = a)
    tilde_D = D
    for (i in 1:m) {
      for (j in 1:m) {
        tilde_D[i, j] = tilde_D[i, j] / (epsilons[i] + epsilons[j])
      }
    }
    valid_scrs = rep(1,50)
    for(k in 2:50){
      valid_scrs[k] = valid_score(k, tilde_D)
    }
    k_selected = which(valid_scrs <= 0.9)[1]
    kmed = cluster::pam(tilde_D, k = k_selected, diss = T)
    CC[id,] = kmed$clustering
  }
  cl_sample = list(KK=apply(CC, 1, max), NN = matrix(0, sample_size, nrow(train)),
                   CC=CC)
  for (j in 1:sample_size) {
    cl_sample$NN[j, 1:cl_sample$KK[j]] = as.numeric(table(CC[j,]))
  }
  m0s = 1:20
  prior_mass = sapply(m0s, function(m0){return(m0 * harmonic_sum(start = m0, end = m0+m-1))})
  prior_mass = abs(prior_mass - mean(cl_sample$KK))
  m0 = which(prior_mass == min(prior_mass))
  save(list = c('cl_sample', 'train', 'test', 'm0', 'test_id'), 
       file = paste('~/CHTC_Downloads/pcba_orig_clust_res_', as.character(idx), '.RData',sep=''))
  
}
rm(dat)


silhouette_score <- function(k){
  tmp_cl = cluster::pam(D, k, diss = T)
  ss <- cluster::silhouette(x = tmp_cl, dist = D)
  mean(ss[, 3])
}
valid_scrs = rep(0,99)
ks = 2:100
for(k in ks){
  valid_scrs[k-1] = valid_score(k, D)
}
plot(ks, type='b', valid_scrs, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
avg_sil <- sapply(ks, silhouette_score)
plot(ks, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
kmed = cluster::pam(D, k = 34, diss = T)
table(kmed$clustering)

adj_mat = matrix(0, nrow(train), nrow(train))
for (i in 1:nrow(train)) {
  for (j in 1:nrow(train)) {
    adj_mat[i,j] = sum(cl_sample$CC[,i]==cl_sample$CC[,j])
  }
}
## visualization
adj_mat = adj_mat / sample_size
ordering = c()
for (k in 1:max(cl_sample$CC[1,])) {
  ordering = c(ordering, which(cl_sample$CC[1,] == k))
}
adj_mat = adj_mat[ordering, ordering]
image(adj_mat, axes = F, xlab = 'Targets', ylab = 'Targets')
mtext(text=seq(0, 100, 10), side=2, line=0.5, at=seq(0, 100, 10) / 101, las=1, cex=1)
mtext(text=seq(0, 100, 10), side=1, line=0.5, at=seq(0, 100, 10) / 101, las=1, cex=1)

### Fast boise
set.seed(817)
sample_size=500
idxs = sample(1:102)#102
## pel1
inform_scores = read.table(paste('~/CHTC_Downloads/PCBA/testid_', as.character(id), '.txt', sep = ''), header = F)
inform = inform_scores$V2[order(inform_scores$V3)[1:1000]]
  
### Overlaps between FDA and PCBA
fda_inform = common_cids
# fda_entropy_compar_result = data.frame('ids' = 1:102, 'roc_fda' = rep(NA, 102),'roc_entropy' = rep(NA, 102), 
#                                        'nef_fda' = rep(NA, 102), 'nef_entropy' = rep(NA, 102))
fda_entropy_compar_result = read.table('~/RAwork/BOISE/BOISE_followup/PCBA_CV/fda_entropy_compar_results.txt', header = T)
fda_entropy_compar_result[, "nef1_entropy"] = rep(NA, 102)
fda_entropy_compar_result[, "nef1_fda"] = rep(NA, 102)
for (id in idxs[1:31]) {
  load(paste('~/CHTC_Downloads/PCBA/pcba_orig_clust_res_', as.character(id), '.RData',sep=''))
  if(sum(test) == 0)
    next
  tmp_inform_fda = which(colnames(train)%in%fda_inform)
  complete_cid = colnames(dat)[which(!is.na(dat[id,]))]
  complete_cid = intersect(complete_cid, fda_inform)
  nA = length(complete_cid) ## valid informer size
  if(nA <= 5){
    print("Few informers")
    next
  }
  # inform_scores = rep(0, ncol(train))
  # for (i in 1:sample_size) {
  #   for (j in 1:ncol(train)) {
  #     x = train[, j]
  #     complete_idx = which(!is.na(x))
  #     prev_entropy = entropy(cl_sample$CC[i, complete_idx])
  #     active_idx = which(x==1)
  #     inactive_idx = which(x==0)
  #     active_rate = mean(x, na.rm = T)
  #     post_entropy = 0
  #     if(length(active_idx) > 0){
  #       post_entropy = post_entropy + active_rate * entropy(cl_sample$CC[i, active_idx])
  #     }
  #     if(length(inactive_idx) > 0){
  #       post_entropy = post_entropy + (1-active_rate) * entropy(cl_sample$CC[i, inactive_idx])
  #     }
  #     mutual_info = prev_entropy - post_entropy
  #     inform_scores[j] = inform_scores[j] + mutual_info
  #   }
  # }
  # tmp_inform_entropy = order(inform_scores, decreasing = T)[1:nA]
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  P = clust_sum(cl_sample,train,sample_size, a, b)
  ## inner function for Evaluate
  Score = rep(0, ncol(train))
  xA = test[tmp_inform_fda]
  nA=length(tmp_inform_fda)
  print(nA)
  m = ncol(train)
  log_post_probs = matrix(0, 1, sample_size)
  post_thetas = matrix(0, sample_size, m)
  for (k in 1:sample_size){
    postls = pel2_beta(P[[k]], x0=train, xA, A=tmp_inform_fda, nT, a, b, m0)
    log_post_probs[k] = postls$log_post_prob
    post_thetas[k, ] = postls$post_theta
  }
  log_post_probs = log_post_probs - max(log_post_probs)
  post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
  Score = post_probs %*% post_thetas
  Score[tmp_inform_fda[which(xA==1)]] = rep(1, sum(xA))
  Score[tmp_inform_fda[which(xA==0)]] = rep(0, nA - sum(xA))
  test = as.vector(test)
  Score = as.vector(Score)
  ## nef
  complete_cid = colnames(dat)[which(!is.na(dat[id,]))]
  nTop = round(length(complete_cid) * 0.01, 0) #nef1
  top = order(Score,decreasing = T)[1:nTop]
  pred_hit = sum(test[top])
  hit = sum(test)
  maxhit = min(hit,nTop)
  nef1 = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
  fda_entropy_compar_result[id, 'nef1_fda'] = nef1
  # nT = as.integer(ncol(train) * 0.1)
  # fda_entropy_compar_result[id, 'roc_fda'] = Evaluate(cl_sample, tmp_inform_fda, 'rocauc',0.1,test,train,
  #                                                    nT, sample_size, a, b, m0)
  # fda_entropy_compar_result[id, 'nef_fda'] = Evaluate(cl_sample, tmp_inform_fda, 'nef',0.1,test,train,
  #                                                    nT, sample_size, a, b, m0)
  # fda_entropy_compar_result[id, 'roc_entropy'] = Evaluate(cl_sample, tmp_inform_entropy, 'rocauc',0.1,test,train,
  #                                                     nT, sample_size, a, b, m0)
  # fda_entropy_compar_result[id, 'nef_entropy'] = Evaluate(cl_sample, tmp_inform_entropy, 'nef',0.1,test,train,
  #                                                     nT, sample_size, a, b, m0)
}
write.table(fda_entropy_compar_result, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/fda_entropy_compar_results.txt',row.names = F)


## entropy
#fast_entropy_result = data.frame('ids' = 1:102, 'roc' = rep(NA, 102), 'nef10' = rep(NA, 102))
fast_entropy_nef_result = read.table('~/RAwork/BOISE/BOISE_followup/PCBA_CV/fast_entropy_nef_1_results.txt', 
                              header = T)
fast_entropy_roc_result = read.table('~/RAwork/BOISE/BOISE_followup/PCBA_CV/fast_entropy_roc_results.txt', 
                                     header = T)
fast_entropy_roc_result[, "roc_500"] = rep(NA, 102)
fast_entropy_nef_result[, "nef_500"] = rep(NA, 102)
nAs = c(100,200,500,1000)
for (id in idxs[1:31]) {
  load(paste('~/CHTC_Downloads/PCBA/pcba_orig_clust_res_', as.character(id), '.RData',sep=''))
  if(sum(test) == 0)
    next
  complete_cid = colnames(dat)[which(!is.na(dat[id,]))]
  inform_scores = rep(0, ncol(train))
  for (i in 1:sample_size) {
    for (j in 1:ncol(train)) {
      x = train[, j]
      complete_idx = which(!is.na(x))
      prev_entropy = entropy(cl_sample$CC[i, complete_idx])
      active_idx = which(x==1)
      inactive_idx = which(x==0)
      active_rate = mean(x, na.rm = T)
      post_entropy = 0
      if(length(active_idx) > 0){
        post_entropy = post_entropy + active_rate * entropy(cl_sample$CC[i, active_idx])
      }
      if(length(inactive_idx) > 0){
        post_entropy = post_entropy + (1-active_rate) * entropy(cl_sample$CC[i, inactive_idx])
      }
      mutual_info = prev_entropy - post_entropy
      inform_scores[j] = inform_scores[j] + mutual_info
    }
  }
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  ## break down Evaluate function to save time!
  P = clust_sum(cl_sample,train,sample_size, a, b)
  for (nA in nAs) {
    inform = order(inform_scores, decreasing = T)[1:nA]
    roc_name = paste('roc_', as.character(nA), sep = '')
    nef_name = paste('nef_', as.character(nA), sep = '')
    fast_entropy_roc_result[id, roc_name] = NA
    fast_entropy_nef_result[id, nef_name] = NA
    ## inner function for Evaluate
    Score = rep(0, ncol(train))
    xA = test[inform]
    nA=length(inform)
    m = ncol(train)
    log_post_probs = matrix(0, 1, sample_size)
    post_thetas = matrix(0, sample_size, m)
    for (k in 1:sample_size){
      postls = pel2_beta(P[[k]], x0=train, xA, A=inform, nT, a, b, m0)
      log_post_probs[k] = postls$log_post_prob
      post_thetas[k, ] = postls$post_theta
    }
    log_post_probs = log_post_probs - max(log_post_probs)
    post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
    Score = post_probs %*% post_thetas
    Score[inform[which(xA==1)]] = rep(1, sum(xA))
    Score[inform[which(xA==0)]] = rep(0, nA - sum(xA))
    test = as.vector(test)
    Score = as.vector(Score)
    ## nef
    nTop = round(length(complete_cid) * 0.01, 0) #nef1
    top = order(Score,decreasing = T)[1:nTop]
    pred_hit = sum(test[top])
    hit = sum(test)
    maxhit = min(hit,nTop)
    nef10 = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
    ## roc
    rocobj = pROC::roc(test,Score)
    rocauc = rocobj$auc
    fast_entropy_roc_result[id, roc_name] = rocauc
    fast_entropy_nef_result[id, nef_name] = nef10
  }
}
write.table(fast_entropy_roc_result, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/fast_entropy_roc_results.txt',row.names = F)
write.table(fast_entropy_nef_result, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/fast_entropy_nef_1_results.txt',row.names = F)

# pel1
fast_pel1_result = read.table('~/RAwork/BOISE/BOISE_followup/PCBA_CV/fast_pel1_results.txt', 
                              header = T)
roc = Evaluate(cl_sample, inform, 'rocauc', 0.1, test, train, nT,sample_size,a,b,m0)
print(roc)
nef10 = Evaluate(cl_sample, inform, 'nef', 0.1, test, train, nT,sample_size,a,b,m0)
print(nef10)
fast_pel1_result$roc[id] = roc
fast_pel1_result$nef10[id] = nef10
write.table(fast_pel1_result, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/fast_pel1_results.txt',row.names = F)

## random selection
rand_roc_result = read.table('~/RAwork/BOISE/BOISE_followup/PCBA_CV/rand_roc_results.txt', header = T)
rand_nef_result = read.table('~/RAwork/BOISE/BOISE_followup/PCBA_CV/rand_nef_results.txt', header = T)
rand_roc_result[, "roc_500"] = rep(NA, 102)
rand_nef_result[, "nef_500"] = rep(NA, 102)
nAs = c(50,100,200,500,1000)
for (id in idxs[1:31]) {
  load(paste('~/CHTC_Downloads/PCBA/pcba_orig_clust_res_', as.character(id), '.RData',sep=''))
  if(sum(test) == 0)
    next
  complete_cid = colnames(dat)[which(!is.na(dat[id,]))]
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  complete_idx = which(!is.na(dat[id,]))
  row_name = colnames(dat)[complete_idx]
  rand_roc_result[id, ] = c(id, rep(0, ncol(rand_roc_result)-1))
  rand_nef_result[id, ] = c(id, rep(0, ncol(rand_nef_result)-1))
  ## break down Evaluate function to save time!
  P = clust_sum(cl_sample,train,sample_size, a, b)
  for(iter in 1:25){
    for (nA in nAs) {
      inform = c()
      while(length(inform) == 0){
        inform = sample(row_name, nA)
        inform = intersect(inform, colnames(train))
      }
      inform = which(colnames(train) %in% inform)
      roc_name = paste('roc_', as.character(nA), sep = '')
      nef_name = paste('nef_', as.character(nA), sep = '')
      ## inner function for Evaluate
      Score = rep(0, ncol(train))
      xA = test[inform]
      nA=length(inform)
      m = ncol(train)
      log_post_probs = matrix(0, 1, sample_size)
      post_thetas = matrix(0, sample_size, m)
      for (k in 1:sample_size){
        postls = pel2_beta(P[[k]], x0=train, xA, A=inform, nT, a, b, m0)
        log_post_probs[k] = postls$log_post_prob
        post_thetas[k, ] = postls$post_theta
      }
      log_post_probs = log_post_probs - max(log_post_probs)
      post_probs = exp(log_post_probs) / (sum(exp(log_post_probs)))
      Score = post_probs %*% post_thetas
      Score[inform[which(xA==1)]] = rep(1, sum(xA))
      Score[inform[which(xA==0)]] = rep(0, nA - sum(xA))
      test = as.vector(test)
      Score = as.vector(Score)
      ## nef
      nTop = round(length(complete_cid) * 0.01, 0)
      top = order(Score,decreasing = T)[1:nTop]
      pred_hit = sum(test[top])
      hit = sum(test)
      maxhit = min(hit,nTop)
      nef10 = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
      ## roc
      rocobj = pROC::roc(test,Score)
      rocauc = rocobj$auc
      rand_roc_result[id, roc_name] = rand_roc_result[id, roc_name] + rocauc
      rand_nef_result[id, nef_name] = rand_nef_result[id, nef_name] + nef10
    }
  }
  rand_roc_result[id, 2:ncol(rand_roc_result)] = rand_roc_result[id, 2:ncol(rand_roc_result)] / 25
  rand_nef_result[id, 2:ncol(rand_nef_result)] = rand_nef_result[id, 2:ncol(rand_nef_result)] / 25
  print(rand_roc_result[id, ])
}
write.table(rand_roc_result, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/rand_roc_results.txt',row.names = F)
write.table(rand_nef_result, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/rand_nef_1_results.txt',row.names = F)


### Plots
library(ggplot2)
idxs = idxs[1:31]
aids = rownames(dat)[idxs]
aids = sapply(aids, function(x){
  tmp = unlist(strsplit(x, 'd'))
  return(tmp[2])
})
fast_entropy_nef_results = read.table('./fast_entropy_nef_1_results.txt', sep = ' ', header=T)
rand_nef_results = read.table('./rand_nef_1_results.txt', sep = ' ', header=T)

legend_title_size = 10
legend_text_size = 10
axis_title_size = 14
axis_text_size = 13
results = data.frame('id'=idxs,'hits' = apply(dat[idxs,],1,function(x) sum(x, na.rm = T)))
results$labels = rep('<100', 31)
results$labels[which(results$hits > 100 & results$hits < 500)] = rep('100~500', 11)
results$labels[which(results$hits > 500)] = rep('>500', 10)
results$labels = as.factor(results$labels)
p1 <- ggplot()+
  geom_point(aes(x = rand_nef_results$nef_100[idxs], y = fast_entropy_nef_results$nef_100[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('random Boise', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of NEF1 for nA=100', color='Hits count')
p2 <- ggplot()+
  geom_point(aes(x = rand_nef_results$nef_200[idxs], y = fast_entropy_nef_results$nef_200[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('random Boise', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of NEF1 for nA=200', color='Hits count')
p3 <- ggplot()+
  geom_point(aes(x = rand_nef_results$nef_500[idxs], y = fast_entropy_nef_results$nef_500[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('random Boise', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of NEF1 for nA=500', color='Hits count')
p4 <- ggplot()+
  geom_point(aes(x = rand_nef_results$nef_1000[idxs], y = fast_entropy_nef_results$nef_1000[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('random Boise', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of NEF1 for nA=1000', color='Hits count')
myGrobs <- list(p1,p2, p3,p4)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 2,ncol = 2)

fast_entropy_results = read.table('./fast_entropy_roc_results.txt', sep = ' ', header=T)
rand_results = read.table('./rand_roc_results.txt', sep = ' ', header=T)
p1 <- ggplot()+
  geom_point(aes(x = rand_results$roc_100[idxs], y = fast_entropy_results$roc_100[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('random Boise', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of ROCAUC for nA=100', color='Hits count')
p2 <- ggplot()+
  geom_point(aes(x = rand_results$roc_200[idxs], y = fast_entropy_results$roc_200[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('random Boise', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of ROCAUC for nA=200', color='Hits count')
p3 <- ggplot()+
  geom_point(aes(x = rand_results$roc_500[idxs], y = fast_entropy_results$roc_500[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('random Boise', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of ROCAUC for nA=500', color='Hits count')
p4 <- ggplot()+
  geom_point(aes(x = rand_results$roc_1000[idxs], y = fast_entropy_results$roc_1000[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('random Boise', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of ROCAUC for nA=1000', color='Hits count')
myGrobs <- list(p1,p2, p3,p4)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 2,ncol = 2)

fda_entropy_compar_result = read.table('./fda_entropy_compar_results.txt', sep = ' ', header=T)
idxs = idxs[1:31]
p1 <- ggplot()+
  geom_point(aes(x = fda_entropy_compar_result$roc_fda[idxs], y = fda_entropy_compar_result$roc_entropy[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('FDA & PCBA intersection', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of ROCAUC for FDA vs. Entropy', color='Hits count' )
p2 <- ggplot()+
  geom_point(aes(x = fda_entropy_compar_result$nef1_fda[idxs], y = fda_entropy_compar_result$nef1_entropy[idxs],color=results$labels))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  scale_x_continuous('FDA & PCBA intersection', limits = c(0.4,1))+
  scale_y_continuous('fast Boise', limits = c(0.4, 1))+
  labs(title = 'Scatterplot of NEF1 for FDA vs. Entropy', color='Hits count')
myGrobs <- list(p1,p2)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 1,ncol = 2)

a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
n = nrow(train)
m0=100
sample_size = 2
### clustering converged even for m0 = 1000!
cl_sample = dpmm_beta(train, a, b, m0= m0, burn_in = 10, sample_size=sample_size, thinning = 5)
save(list = c('cl_sample', 'm0', 'sample_size'), file = 'pcba_clust_reduced.RData')
load('pcba_clust_reduced.RData')
sample_size = 2
adj_mat = matrix(0, nrow(train), nrow(train))
for (i in 1:nrow(train)) {
  for (j in 1:nrow(train)) {
    adj_mat[i,j] = sum(cl_sample$CC[,i]==cl_sample$CC[,j])
  }
}
adj_mat = adj_mat / sample_size
ordering = c()
for (k in 1:cl_sample$KK[1]) {
  ordering = c(ordering, which(cl_sample$CC[1,] == k))
}
adj_mat = adj_mat[ordering, ordering]
image(adj_mat, axes = F, xlab = 'Targets', ylab = 'Targets')
mtext(text=seq(0, 100, 10), side=2, line=0.5, at=seq(0, 100, 10) / 101, las=1, cex=1)
mtext(text=seq(0, 100, 10), side=1, line=0.5, at=seq(0, 100, 10) / 101, las=1, cex=1)

acf(cl_sample$KK)
which(cl_sample$CC[1,] != cl_sample$CC[2,])
cl_sample$NN[1,1:26]

grp = which(cl_sample$CC[1,]==1)
tmp = dat[grp,]
tmp_means = colMeans(tmp, na.rm = T)
sort(tmp_means,decreasing = T)[10]
tmp[, order(tmp_means,decreasing = T)[1:20]]

### analyze the convergence of clustering
cl = Initial_beta(train, m0)
count = 0
cl_next = Update_beta(cl, train, a, b, m0)
while (count < 2000 & sum(cl$C != cl_next$C) > 0) {
  cl = cl_next
  cl_next = Update_beta(cl, train, a, b, m0)
  count = count + 1
}
print(count)
