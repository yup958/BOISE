setwd('~/RAwork/BOISE/BOISE_followup/')
load('fda_data_rearranged.RData')
block_res$row_block = block_res$row_block[1:5,]
block_res$col_block = block_res$col_block[1:5,]
active_ind = apply(block_res$col_block, 2, sum)
active_ind = which(active_ind > 0)
cpds = colnames(block_res$col_block[, active_ind])
dat = dat[, cpds] ### Missing rate on residual is 92.9%
counts = apply(dat, 1, function(x){return(sum(!is.na(x)))})
sort(counts, decreasing = T)[50] ## max complete target: 826/933; top50: 669, top100: 622

### For max complete data, with 688 x 508, penalization 20.
# load('max_complete_data.RData')
# A = dat[, setdiff(colnames(dat),largest_block$CID)] 
## 688x508; Missing rate on residual is 87.6%. Top 100 complete targets has at least 507 complete compounds!
# dat = dat[,largest_block$CID]

## for 688x508 data, 35~45 clusters, 
## with top 6 clusters of size >50; largest cluster 120~150; 24~28 clusters of size < 10; 11~14 clusters <= 2

library(BOISE)
set.seed(2)
a = rep(mean(dat, na.rm = T), ncol(dat))
b = 1 - a

## helper function to find m0
harmonic_sum <- function(start, end){
  res = 0
  for(k in start:end){
    res = res + 1/k
  }
  return(res)
}

m0_Find <- function(x0, lower, upper){
  a = rep(mean(x0), ncol(x0))
  b = 1 - a
  n = nrow(x0)
  min_diff = 10000
  best_m0 = 1
  
  for(m0 in lower:upper){
    prior_m0_expect = m0 * harmonic_sum(start = m0, end = m0+n-1)
    cl_sample = dpmm_beta(x0, a, b, m0, burn_in = 1000, sample_size=100, thinning = 10)
    diff = abs(mean(cl_sample$KK) - prior_m0_expect)
    if(diff < min_diff){
      min_diff = diff
      best_m0 = m0
    }
  }
  return(best_m0)
}

m0 = m0_Find(dat, lower = 1, upper = 15) 
## Clustering whole matrix (688 x 933): select m0 = 9.
## Clustering each group of compounds (1-5):
## Grp1 m0 = 13; Grp2 m0 = 10; Grp3 m0 = 8; Grp4 m0 = 3; Grp5 m0 = 2.
m0 = 9

### BOISE on the whole data set
sample_size = 100
system.time({
  cl_sample = dpmm_beta(dat, a, b, m0= m0, burn_in = 1000, sample_size=sample_size, thinning = 10)
}) ## 5h15min
save(list = c('block_res', 'dat', 'cl_sample'), file = 'clust_res_whole_mat.RData')
## for 688 x 933 data, clustering in whole, m0=9; end up with 36-46 clusters
counts = apply(cl_sample$NN, 2, mean)
counts = counts[which(counts>1)] # 35 common clusters. 12 clusters > 10; 23 clusters < 9, 15 clusters < 3.
print(counts)



### BOISE on separate groups of 688 x 933 data:
sample_size = 100
m0_selections = c(13, 10, 8, 3, 2)
cl_samples = list(1,2,3,4,5)
for(k in 1:5){
  active_cid = which(block_res$col_block[k, ] > 0)
  cpds = colnames(block_res$col_block[, active_cid])
  A = dat[, cpds] 
  a = rep(mean(A, na.rm = T), ncol(A))
  b = 1 - a
  m0 = m0_selections[k]
  cl_samples[[k]] = dpmm_beta(A, a, b, m0, burn_in = 1000, sample_size=sample_size, thinning = 10)
}

### Posterior predictive check on non-missing entries?
load('clust_res_whole_mat.RData')
generate_posterior_parameters <- function(cl, x0, alpha, beta){
  m = nrow(x0)
  n = ncol(x0)
  K = cl$K
  Phi = matrix(0, K, n)
  for (ck in 1:K) {
    target = which(cl$C== ck)
    if(length(target)==1){
      missing = which(is.na(x0[target, ]))
      ak = alpha + x0[target, ]
      bk = beta + (1 - x0[target,])
      ak[missing] = alpha[missing]
      bk[missing] = beta[missing]
    }else{
      ak = alpha + colSums(x0[target,], na.rm = T)
      bk = beta + colSums(1 - x0[target,], na.rm = T)
    }
    for (j in 1:n) {
      Phi[ck, j] = rbeta(1, ak[j], bk[j])
    }
  }
  return(Phi)
}
generate_fake_data <- function(cl, Phi){
  n = ncol(Phi)
  m = length(cl$C)
  Y = matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    ck = cl$C[i]
    P = Phi[ck,]
    Y[i, ] = rbinom(n, 1, P)
  }
  return(Y)
}
log_likelihood <- function(cl, Phi, mat, x0){
  m = nrow(x0)
  neg_log_lik = 0
  for(i in 1:m){
    ck = cl$C[i]
    P = Phi[ck, ]
    non_missing_idx = which(!is.na(x0[i, ]))
    neg_log_lik = neg_log_lik - sum(dbinom(mat[i, non_missing_idx], 1, P[non_missing_idx], log = TRUE))
  }
  return(neg_log_lik)
}
## for clustering on whole matrix
sample_size=100
iter = 100
count = 0
a = rep(mean(dat, na.rm = T), ncol(dat))
b = 1 - a
# post pred check
for (i in 1:sample_size) {
  cl = list(K = cl_sample$KK[i], N = cl_sample$NN[i, ], C = cl_sample$CC[i, ])
  for(j in 1:iter){
    post_phi = generate_posterior_parameters(cl, dat, a, b)
    Y = generate_fake_data(cl, post_phi)
    real_log_lik = log_likelihood(cl, post_phi, mat = dat, x0 = dat)
    fake_log_lik = log_likelihood(cl, post_phi, mat = Y, x0 = dat)
    if(fake_log_lik > real_log_lik){
      count = count + 1
    }
  }
}
print(count / sample_size / iter) ## 0.0028 for clustering on whole matrix, B = 10000.
# average log likelihood
real_log_lik = 0
for (i in 1:sample_size) {
  cl = list(K = cl_sample$KK[i], N = cl_sample$NN[i, ], C = cl_sample$CC[i, ])
  for(j in 1:iter){
    post_phi = generate_posterior_parameters(cl, dat, a, b)
    real_log_lik = real_log_lik + log_likelihood(cl, post_phi, mat = dat, x0 = dat)
  }
}
print(real_log_lik / sample_size / iter) ## 21037.1


## for clustering in groups
load('clust_res_in_grps.RData')
sample_size=100
iter = 100
count = 0
# post pred check
for (i in 1:sample_size) {
  for(j in 1:iter){
    real_log_lik = 0
    fake_log_lik = 0
    
    for(k in 1:5){
      cl = list(K = cl_samples[[k]]$KK[i], N = cl_samples[[k]]$NN[i, ], C = cl_samples[[k]]$CC[i, ])
      active_cid = which(block_res$col_block[k, ] > 0)
      cpds = colnames(block_res$col_block[, active_cid])
      A = dat[, cpds] 
      a = rep(mean(dat, na.rm = T), ncol(A))
      b = 1 - a
      post_phi = generate_posterior_parameters(cl, A, a, b)
      Y = generate_fake_data(cl, post_phi)
      real_log_lik = real_log_lik + log_likelihood(cl, post_phi, mat = A, x0 = A)
      fake_log_lik = fake_log_lik + log_likelihood(cl, post_phi, mat = Y, x0 = A)
    }
    
    if(fake_log_lik > real_log_lik){
      count = count + 1
    }
  }
}
print(count / sample_size / iter) ## 0.1317 for clustering on whole matrix, B = 10000.
# avg log likelihood
real_log_lik = 0
for (i in 1:sample_size) {
  for(j in 1:iter){
    for(k in 1:5){
      cl = list(K = cl_samples[[k]]$KK[i], N = cl_samples[[k]]$NN[i, ], C = cl_samples[[k]]$CC[i, ])
      active_cid = which(block_res$col_block[k, ] > 0)
      cpds = colnames(block_res$col_block[, active_cid])
      A = dat[, cpds] 
      a = rep(mean(dat, na.rm = T), ncol(A))
      b = 1 - a
      post_phi = generate_posterior_parameters(cl, A, a, b)
      real_log_lik = real_log_lik + log_likelihood(cl, post_phi, mat = A, x0 = A)
    }
  }
}
print(real_log_lik / sample_size / iter) ## 18395.15


counts = apply(cl_samples[[2]]$NN,2,mean)
counts = counts[counts>1]
counts
