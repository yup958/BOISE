setwd("~/RAwork/BOISE/BOISE_followup/PCBA_prospective/")
test = read.csv('./training_df_single_fold.csv.gz')
test = test[,c('rdkit.SMILES', 'PriA.SSB.AS.Activity')]
train = read.csv('~/Downloads/condensed_subset_PCBA128.csv.gz')
train = train %>% distinct(rdkit.SMILES, .keep_all = TRUE)
## inner join on SMILES
dat = merge(x=train, y = test, by = 'rdkit.SMILES', all = FALSE)
rownames(dat) = dat[,'PUBCHEM_CID']
dat = dat[, !colnames(dat) %in% c("rdkit.SMILES","Morgan.FP_2_1024","Cluster.ID", "PUBCHEM_CID")]
train = dat[,1:102]
test = dat[,103]
names(test) = rownames(dat)
train = t(train)
print(sum(is.na(test))) ### no missing values on PriA_SSB

counts = colSums(train, na.rm = T)
train = train[, which(counts > 0)]

# function to compute validity score
valid_score <- function(k, D=D){
  kmed = cluster::pam(D, k, diss = T)
  # among cluster distance
  inter = 0
  for(i in 1:(k-1)){
    medoid1 = kmed$medoids[i]
    for (j in (i+1):k) {
      medoid2 = kmed$medoids[j]
      inter = inter + (D[medoid1, medoid2])
    }
  }
  inter = inter / (k*(k-1) / 2)
  # within cluster distance
  m = nrow(D)
  intra = rep(0, m)
  for (i in 1:m) {
    label = kmed$clustering[i]
    med = kmed$medoids[label]
    intra[i] = (D[i, med])
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
sample_size = 500

### Approximate clustering with random distance using Hamming / Manhatton dist or Jaccard dist
m = nrow(train)
D = matrix(0,m,m)
for (i in 1:m) {
  for (j in 1:m) {
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
  m = 102
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
optim_obj = nlminb(start = c(10,1000), objective = log_obj, lower = 0, upper = 1e6,control = list(eval.max = 600))
a = sum(optim_obj$par) 

CC = matrix(0, sample_size, nrow(train))
set.seed(817)
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
  k_selected = which(valid_scrs <= 0.7)[1]
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
#save(list = c('cl_sample', 'train', 'test', 'm0', 'test_id'), 
#     file = paste('~/CHTC_Downloads/pcba_orig_clust_res_', as.character(idx), '.RData',sep=''))

### Fast boise
set.seed(817)
sample_size=500
library(entropy)
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
nA=1000
inform = order(inform_scores, decreasing = T)[1:nA]
## inner function for Evaluate
Score = rep(0, ncol(train))
truncated_test = test[which(names(test)%in%colnames(train))]
xA = truncated_test[inform]
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
names(Score) = colnames(train)
pred_scores = rep(min(Score), length(test))
names(pred_scores) = names(test)
for (i in 1:length(Score)) {
  pred_scores[which(names(pred_scores) == names(Score)[i])] = Score[i]
}
## nef
nTop = round(length(test) * 0.01, 0) #nef1
top = order(pred_scores,decreasing = T)[1:nTop]
pred_hit = sum(test[top])
hit = sum(test)
maxhit = min(hit,nTop)
nef1 = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
print(nef1)
## roc
rocobj = pROC::roc(test,pred_scores)
rocauc = rocobj$auc
print(rocauc)
fast_boise_res = c(nef1, rocauc) ## for nA = 1000, nef1 = 0.622, rocauc = 0.890

### rand boise
## break down Evaluate function to save time!
a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
P = clust_sum(cl_sample,train,sample_size, a, b)
nA=253
pubids = names(test)
set.seed(817)
m = ncol(train)
rand_roc = 0
rand_nef1 = 0
for(iter in 1:25){
  inform = c()
  while(length(inform) == 0){
    inform = sample(pubids, nA)
    inform = intersect(inform, colnames(train))
  }
  inform = which(colnames(train) %in% inform)
  ## inner function for Evaluate
  Score = rep(0, ncol(train))
  truncated_test = test[which(names(test)%in%colnames(train))]
  xA = truncated_test[inform]
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
  names(Score) = colnames(train)
  pred_scores = rep(min(Score), length(test))
  names(pred_scores) = names(test)
  for (i in 1:length(Score)) {
    pred_scores[which(names(pred_scores) == names(Score)[i])] = Score[i]
  }
  ## nef
  nTop = round(length(test) * 0.01, 0) #nef1
  top = order(pred_scores,decreasing = T)[1:nTop]
  pred_hit = sum(test[top])
  hit = sum(test)
  maxhit = min(hit,nTop)
  nef1 = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2
  ## roc
  rocobj = pROC::roc(test,pred_scores)
  rocauc = rocobj$auc
  rand_roc = rand_roc + rocauc
  rand_nef1 = rand_nef1 + nef1
}
rand_roc = rand_roc / 25
rand_nef1 = rand_nef1 / 25
rand_boise_res = c(rand_nef1, rand_roc) ## for nA=1000, nef1 = 0.533, rocauc = 0.893

### overlaps with FDA
load('~/RAwork/BOISE/BOISE_followup/FDA_CV/fda_data.RData')
fda_cids = colnames(dat)
fda_cids = sapply(fda_cids, function(s){return(unlist(strsplit(s, '_'))[2])})
common_cids = intersect(fda_cids, names(test)) ##253 common compounds in both FDA and LC
# For fast boise with nA=253
nA=253
inform = order(inform_scores, decreasing = T)[1:nA] ## nef1 0.553, rocauc 0.862
# For FDA inform
inform = which(colnames(train) %in% common_cids) ## nef1 0.536, rocauc 0.838

## visualization
a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
P = clust_sum(cl_sample,train,sample_size, a, b)
nA=253
inform = order(inform_scores, decreasing = T)[1:nA]
truncated_test = test[which(names(test)%in%colnames(train))]
xA = truncated_test[inform]
post_adj_mat = matrix(0, sample_size, nrow(train) + 1)
for (i in 1:sample_size) {
  lp = apply(P[[i]][,1:ncol(train)], 1, function(q){
    return(sum(xA * log(q[inform])) + sum((1-xA) * log(1 - q[inform])))
  })
  lnew = sum(xA * log(a[inform] / (a[inform] + b[inform]))) + 
    sum((1 - xA) * log(b[inform] / (a[inform] + b[inform])))
  lp = c(lp, lnew) # P(xA|i in ck, C, x0)
  
  prob_masses = c(P[[i]][,ncol(train)+1], m0)
  log_prob_prior = log(prob_masses / sum(prob_masses)) #P(i in ck|C, x0)
  
  log_prob_prior = log_prob_prior + lp
  c = max(log_prob_prior)
  log_prob_prior = log_prob_prior - c
  log_post_probs = c + log(sum(exp(log_prob_prior))) #P(xA|C, x0)
  
  weights = lp + log_prob_prior - log_post_probs
  c = max(weights)
  weights = weights - c
  weights = exp(weights) / sum(exp(weights))
  post_adj_mat[i, nrow(train) + 1] = weights[length(weights)]
  post_adj_mat[i, 1:nrow(train)] = sapply(cl_sample$CC[i,], function(c) return(weights[c]))
}

adj_mat = matrix(0, nrow(train), nrow(train))
for (i in 1:nrow(train)) {
  for (j in 1:nrow(train)) {
    adj_mat[i,j] = sum(cl_sample$CC[,i]==cl_sample$CC[,j])
  }
}
adj_mat = adj_mat / sample_size
## add a "fake" target, for unknown new table in CRP
adj_mat = cbind(adj_mat, rep(0, 102))
adj_mat = rbind(adj_mat, c(rep(0, 102), 1))

wt_adj_mat = matrix(0, nrow(train)+1, nrow(train)+1)
for (k in 1:sample_size) {
  for (i in 1:nrow(train)) {
    for (j in 1:nrow(train)) {
      wt_adj_mat[i,j] = wt_adj_mat[i,j] + 
        (cl_sample$CC[k,i]==cl_sample$CC[k,j]) * post_adj_mat[k, i] ##or [k,j]
    }
  }
  wt_adj_mat[nrow(train)+1, nrow(train)+1] = 
    wt_adj_mat[nrow(train)+1, nrow(train)+1] + post_adj_mat[k, nrow(train)+1]
}
wt_adj_mat = wt_adj_mat / max(wt_adj_mat)

## if we want to plot similarity of the novel target to existing clusters:
# post_adj_mat = colMeans(post_adj_mat)
# post_adj_mat = as.matrix(post_adj_mat)
# wt_mat = post_adj_mat %*% t(post_adj_mat)
# wt_adj_mat = adj_mat * wt_mat
# wt_adj_mat = wt_adj_mat / max(wt_adj_mat)
 
ordering = c()
for (k in 1:max(cl_sample$CC[1,])) {
  ordering = c(ordering, which(cl_sample$CC[1,] == k))
}
ordering = c(ordering, 103)
adj_mat = adj_mat[ordering, ordering]
par(mfrow = c(1,2))
image(adj_mat, axes = F, xlab = 'Targets', ylab = 'Targets')
mtext(text=seq(0, 100, 10), side=2, line=0.5, at=seq(0, 100, 10) / 103, las=1, cex=1)
mtext(text=seq(0, 100, 10), side=1, line=0.5, at=seq(0, 100, 10) / 103, las=1, cex=1)

wt_adj_mat = wt_adj_mat[ordering, ordering]
image(wt_adj_mat, axes = F, xlab = 'Targets', ylab = 'Targets')
mtext(text=seq(0, 100, 10), side=2, line=0.5, at=seq(0, 100, 10) / 103, las=1, cex=1)
mtext(text=seq(0, 100, 10), side=1, line=0.5, at=seq(0, 100, 10) / 103, las=1, cex=1)

## unknown target (103) is most similar to PriA-SSB
print(order(post_adj_mat, decreasing = T))
most_sim_tgts = ordering[which(wt_adj_mat[79,] > 0.1)] ##103, 31, 77, 69, 17
most_sim_tgts = 77
most_sim_sub_matrix = t(dat)[most_sim_tgts,]
# compare similarity between unknown target and test
avg_active = rowMeans(dat, na.rm = T)
print(sum(test * avg_active)) #3.381, seems good
rocauc = pROC::roc(test,avg_active)$auc
print(rocauc) #0.9546 rocauc!
nTop = round(length(test) * 0.01, 0) #nef1
top = order(avg_active,decreasing = T)[1:nTop]
pred_hit = sum(test[top])
hit = sum(test)
maxhit = min(hit,nTop)
nef1 = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2 #0.6565

# compare similarity between most similar cluster and test
#avg_active = colMeans(most_sim_sub_matrix, na.rm = T)
avg_active = most_sim_sub_matrix
avg_active[which(is.na(avg_active))] = rowMeans(dat[which(is.na(avg_active)),], na.rm = T)
print(sum(test * avg_active))
rocauc = pROC::roc(test,avg_active)$auc
print(rocauc)
nTop = round(length(test) * 0.01, 0) #nef1
top = order(avg_active,decreasing = T)[1:nTop]
pred_hit = sum(test[top])
hit = sum(test)
maxhit = min(hit,nTop)
nef1 = ((pred_hit/nTop - hit/ncol(train)) / (maxhit/nTop - hit/ncol(train)) + 1)/2 
print(nef1)
## for tgt id 31, nef1 0.562; for tgt id 77, 0.519; for tgt id 69, 0.511...; for id 17, 0.493
