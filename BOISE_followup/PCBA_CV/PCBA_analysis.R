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

idxs = sample(1:102)#102
sample_size = 100
for(idx in idxs[11:30]){
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
  optim_obj = nlminb(start = c(300,300), objective = log_obj, lower = 0, upper = 1e6,control = list(eval.max = 600))
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
    k_selected = which(valid_scrs <= 0.8)[1]
    kmed = cluster::pam(tilde_D, k = k_selected, diss = T)
    CC[id,] = kmed$clustering
  }
  cl_sample = list(KK=apply(CC, 1, max), NN = matrix(0, sample_size, nrow(train)),
                   CC=CC)
  for (j in 1:sample_size) {
    cl_sample$NN[j, 1:cl_sample$KK[j]] = as.numeric(table(CC[j,]))
  }
  m0 = 2
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
    adj_mat[i,j] = sum(CC[,i]==CC[,j])
  }
}
## visualization
adj_mat = adj_mat / sample_size
ordering = c()
for (k in 1:max(CC[1,])) {
  ordering = c(ordering, which(CC[1,] == k))
}
adj_mat = adj_mat[ordering, ordering]
image(adj_mat, axes = F, xlab = 'Targets', ylab = 'Targets')
mtext(text=seq(0, 100, 10), side=2, line=0.5, at=seq(0, 100, 10) / 101, las=1, cex=1)
mtext(text=seq(0, 100, 10), side=1, line=0.5, at=seq(0, 100, 10) / 101, las=1, cex=1)

### Fast boise
set.seed(817)
sample_size=100
idxs = sample(1:102)#102
id = idxs[15]
load(paste('~/CHTC_Downloads/PCBA/pcba_orig_clust_res_', as.character(id), '.RData',sep=''))
a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
P = clust_sum(cl_sample, train, sample_size, a, b)
m = nrow(train)
n = ncol(train)
nA = 1
nT = as.integer(ncol(train) * 0.1)
## pel1
inform_scores = read.table(paste('~/CHTC_Downloads/PCBA/testid_', as.character(id), '.txt', sep = ''), header = F)
inform = inform_scores$V2[order(inform_scores$V3)[1:1000]]
  

## entropy
#fast_entropy_result = data.frame('ids' = 1:102, 'roc' = rep(NA, 102), 'nef10' = rep(NA, 102))
fast_entropy_result = read.table('~/RAwork/BOISE/BOISE_followup/PCBA_CV/fast_entropy_results.txt', 
                              header = T)
nA = 200
for (id in idxs[1:15]) {
  load(paste('~/CHTC_Downloads/PCBA/pcba_orig_clust_res_', as.character(id), '.RData',sep=''))
  if(sum(test) == 0)
    next
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
  inform = order(inform_scores, decreasing = T)[1:nA]
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  roc = Evaluate(cl_sample, inform, 'rocauc', 0.1, test, train, nT,sample_size,a,b,m0)
  print(roc)
  nef10 = Evaluate(cl_sample, inform, 'nef', 0.1, test, train, nT,sample_size,a,b,m0)
  print(nef10)
  fast_entropy_result$roc[id] = roc
  fast_entropy_result$nef10[id] = nef10
}
write.table(fast_roc_results, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/fast_entropy_roc_results.txt',row.names = F)
write.table(fast_nef_results, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/fast_entropy_nef_results.txt',row.names = F)

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
rand_result = read.table('~/RAwork/BOISE/BOISE_followup/PCBA_CV/rand_results.txt', header = T)
for (id in idxs[1:15]) {
  load(paste('~/CHTC_Downloads/PCBA/pcba_orig_clust_res_', as.character(id), '.RData',sep=''))
  if(sum(test) == 0)
    next
  a = rep(mean(train, na.rm = T), ncol(train))
  b = 1 - a
  roc = 0
  nef10 = 0
  for(iter in 1:25){
    inform = sample(1:ncol(train), nA)
    roc = roc + Evaluate(cl_sample, inform, 'rocauc', 0.1, test, train, nT,sample_size,a,b,m0)
    nef10 = nef10 + Evaluate(cl_sample, inform, 'nef', 0.1, test, train, nT,sample_size,a,b,m0)
  }
  roc = roc / 25
  nef10 = nef10 / 25
  print(roc)
  print(nef10)
  rand_result$roc[id] = roc
  rand_result$nef10[id] = nef10
}
write.table(rand_roc_results, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/rand_roc_results.txt',row.names = F)
write.table(rand_nef_results, file = '~/RAwork/BOISE/BOISE_followup/PCBA_CV/rand_nef_results.txt',row.names = F)


### Plots
library(ggplot2)
idxs = idxs[1:15]
aids = rownames(dat)[idxs]
aids = sapply(aids, function(x){
  tmp = unlist(strsplit(x, 'd'))
  return(tmp[2])
})
fast_entropy_results = read.table('./fast_entropy_nef_results.txt', sep = ' ', header=T)
rand_results = read.table('./rand_nef_results.txt', sep = ' ', header=T)

legend_title_size = 10
legend_text_size = 10
axis_title_size = 14
axis_text_size = 13
p1 <- ggplot()+
  geom_point(aes(x = rand_results$nef_50[idxs], y = fast_entropy_results$nef_50[idxs]))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of NEF10 for nA=50', x = 'random Boise', y = 'fast Boise')
p2 <- ggplot()+
  geom_point(aes(x = rand_results$nef_100[idxs], y = fast_entropy_results$nef_100[idxs]))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of NEF10 for nA=100', x = 'random Boise', y = 'fast Boise')
p3 <- ggplot()+
  geom_point(aes(x = rand_results$nef_200[idxs], y = fast_entropy_results$nef_200[idxs]))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of NEF10 for nA=200', x = 'random Boise', y = 'fast Boise')
p4 <- ggplot()+
  geom_point(aes(x = rand_results$nef_1000[idxs], y = fast_entropy_results$nef_1000[idxs]))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of NEF10 for nA=1000', x = 'random Boise', y = 'fast Boise')
myGrobs <- list(p1,p2, p3,p4)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 2,ncol = 2)

fast_entropy_results = read.table('./fast_entropy_roc_results.txt', sep = ' ', header=T)
rand_results = read.table('./rand_roc_results.txt', sep = ' ', header=T)
p1 <- ggplot()+
  geom_point(aes(x = rand_results$roc_50[idxs], y = fast_entropy_results$roc_50[idxs]))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of ROCAUC for nA=50', x = 'random Boise', y = 'fast Boise')
p2 <- ggplot()+
  geom_point(aes(x = rand_results$roc_100[idxs], y = fast_entropy_results$roc_100[idxs]))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of ROCAUC for nA=100', x = 'random Boise', y = 'fast Boise')
p3 <- ggplot()+
  geom_point(aes(x = rand_results$roc_200[idxs], y = fast_entropy_results$roc_200[idxs]))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of ROCAUC for nA=200', x = 'random Boise', y = 'fast Boise')
p4 <- ggplot()+
  geom_point(aes(x = rand_results$roc_1000[idxs], y = fast_entropy_results$roc_1000[idxs]))+
  geom_abline(slope = 1, intercept = 0, color = 'red')+
  labs(title = 'Scatterplot of ROCAUC for nA=1000', x = 'random Boise', y = 'fast Boise')
myGrobs <- list(p1,p2, p3,p4)
gridExtra::grid.arrange(grobs = myGrobs, nrow = 2,ncol = 2)

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
