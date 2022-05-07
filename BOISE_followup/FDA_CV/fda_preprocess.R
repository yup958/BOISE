setwd('~/RAwork/BOISE/BOISE_followup/FDA_CV/')
set.seed(817)
load('fda_data.RData')
# cpds = colnames(dat)
# write.table(cpds, file = 'Compound_IDS.txt',col.names = F, row.names = F)

### Sample test ids
test_id = c()
while (length(test_id) < 60) {
  id = sample(rownames(dat),1)
  len_no_miss = sum(!is.na(dat[id,]))
  if(len_no_miss > 400){
    test_id = c(test_id, id)
  }
}
apply(dat[test_id,],1, function(x){return(sum(!is.na(x)))})
apply(dat[test_id,],1, function(x){return(sum(x, na.rm=T))})
write.table(test_id, file = 'Test_IDS.txt',col.names = F, row.names = F)

### read distance matrix, use chemical clustering
dist = data.table::fread('Similarity_mat_FDA.csv.gz')
rownames(dist) = dist$V1
dist$V1 = NULL
dist = 1 - dist
cl = cluster::diana(x = dist, diss = T, stop.at.k = 20, keep.diss = T)

# hcd = as.dendrogram(cl)
# tmp = cut(hcd, h = 0.95)
# plot(tmp$upper)
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

## for original boise
set.seed(817)
sample_size=100
iter = 1
count = 0
id=2
load(paste('~/CHTC_Downloads/FDA_cv/orig_clust_res_', as.character(id), '.RData',sep=''))
a = rep(mean(train, na.rm = T), ncol(train))
b = 1 - a
# post pred check
for (i in 1:sample_size) {
  cl = list(K = cl_sample$KK[i], N = cl_sample$NN[i, ], C = cl_sample$CC[i, ])
  for(j in 1:iter){
    post_phi = generate_posterior_parameters(cl, train, a, b)
    Y = generate_fake_data(cl, post_phi)
    #Y = matrix(0, nrow=nrow(train), ncol = ncol(train))
    real_log_lik = log_likelihood(cl, post_phi, mat = train, x0 = train)
    fake_log_lik = log_likelihood(cl, post_phi, mat = Y, x0 = train)
    if(fake_log_lik > real_log_lik){
      count = count + 1
    }
  }
}
print(count / sample_size / iter) 

## for clustering in groups
set.seed(817)
sample_size=100
iter = 1
id=20
load(paste('~/CHTC_Downloads/FDA_cv/testid_', as.character(id), '_block.RData',sep=''))
count = 0
# post pred check
for (i in 1:sample_size) {
  for(j in 1:iter){
    real_log_lik = 0
    fake_log_lik = 0
    col_clust_num = length(block)
    for(k in 1:col_clust_num){
      cl_sample = list(KK=apply(block[[k]], 1, max), NN = matrix(0, sample_size, ncol(block[[k]])),
                       CC=block[[k]])
      for (ll in 1:sample_size) {
        cl_sample$NN[ll, 1:cl_sample$KK[ll]] = as.numeric(table(block[[k]][ll,]))
      }
      tmp_cl = list(K = cl_sample$KK[i], N = cl_sample$NN[i, ], C = cl_sample$CC[i, ])
      cpds = cl$cid[which(cl$C==k)]
      A = train[, cpds] 
      a = rep(mean(A, na.rm = T), ncol(A))
      b = 1 - a
      post_phi = generate_posterior_parameters(tmp_cl, A, a, b)
      Y = generate_fake_data(tmp_cl, post_phi)
      real_log_lik = real_log_lik + log_likelihood(tmp_cl, post_phi, mat = A, x0 = A)
      fake_log_lik = fake_log_lik + log_likelihood(tmp_cl, post_phi, mat = Y, x0 = A)
    }
    
    if(fake_log_lik > real_log_lik){
      count = count + 1
    }
  }
}
print(count / sample_size / iter) 

