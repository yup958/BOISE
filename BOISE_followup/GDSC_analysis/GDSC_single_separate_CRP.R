setwd('~/RAwork/BOISE/BOISE_followup/')
set.seed(2)
load('GDSC.RData')
#load('GDSC.RData')

### Single CRP Posterior Predictive check
m0 = 30 ## select 30 for GDSC
sample_size = 200
a = rep(mean(GDSC, na.rm = T), ncol(GDSC))
b = 1 - a
system.time({
  cl_sample = dpmm_beta(GDSC, a, b, m0= m0, burn_in = 1500, sample_size=sample_size, thinning = 15)
}) 
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
iter = 100
count = 0

# post pred check
for (i in 1:sample_size) {
  cl = list(K = cl_sample$KK[i], N = cl_sample$NN[i, ], C = cl_sample$CC[i, ])
  for(j in 1:iter){
    post_phi = generate_posterior_parameters(cl, GDSC, a, b)
    Y = generate_fake_data(cl, post_phi)
    real_log_lik = log_likelihood(cl, post_phi, mat = GDSC, x0 = GDSC)
    fake_log_lik = log_likelihood(cl, post_phi, mat = Y, x0 = GDSC)
    if(fake_log_lik > real_log_lik){
      count = count + 1
    }
  }
}
print(count / sample_size / iter) ## 1 / 20000 for single CRP B=20000
# average log likelihood
real_log_lik = 0
for (i in 1:sample_size) {
  cl = list(K = cl_sample$KK[i], N = cl_sample$NN[i, ], C = cl_sample$CC[i, ])
  for(j in 1:iter){
    post_phi = generate_posterior_parameters(cl, GDSC, a, b)
    real_log_lik = real_log_lik + log_likelihood(cl, post_phi, mat = GDSC, x0 = GDSC)
  }
}
print(real_log_lik / sample_size / iter) ## 8571.25



### Separate CRPs

# Spectral clustering
jaccard_dist <- function(x, y){
  positive_x = which(x == 1)
  positive_y = which(y == 1)
  inter = intersect(positive_x, positive_y)
  uni = union(positive_x, positive_y)
  if(length(uni) == 0){
    return(0)
  } else{
    return(1 -  length(inter)/length(uni))
  }
}

hamming_dist <- function(x, y){
  ## x and y are both binary vectors
  return(sum(abs(x-y)))
}
"%^%" <- function(M, power)
  with(eigen(M), vectors %*% (values^power * solve(vectors)))
spectral_clust <- function(x0, k){
  n = ncol(x0)
  dist_mat = matrix(0, n, n)
  for (j in 1:n) {
    y = x0[ , j]
    # dist_mat[ , j] = apply(x0, 2, function(x){return(jaccard_dist(x, y))})
    dist_mat[ , j] = apply(x0, 2, function(x){return(hamming_dist(x, y))})
  }
  
  #A = 1 - dist_mat    #Affinity Matrix for Jaccard
  A = 1 - dist_mat / n
  D = diag(apply(A, 1, sum))
  #L = D - A    ##Unnormalized Laplacian
  L = diag(rep(1, n)) - solve(D) %*% A    ##Simple Laplacian
  #L <- diag(rep(1, n)) - (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))  # Normalized Laplacian
  eig <- eigen(L, symmetric = TRUE)
  Z <- eig$vectors[,(ncol(eig$vectors)-k + 1):ncol(eig$vectors)]
  km <- kmeans(Z, iter.max = 200, centers = k, nstart = 20)
  return(km)
}

for (k in 2:10) {
  print(k)
  cl = spectral_clust(x0 = GDSC, k = k)
  print(cl$tot.withinss / cl$betweenss) ## Let's select k = 6
}

cl = spectral_clust(x0 = GDSC, k = 2)
col_block = matrix(0,nrow = max(cl$cluster), ncol = length(cl$cluster))
colnames(col_block) = colnames(GDSC)
for (j in 1:ncol(col_block)) {
  col_block[cl$cluster[j], j] = 1
}
save(list = c('col_block', 'GDSC'), file = 'GDSC_hamming_dist.RData')
# Jaccard dist: cluster sizes are c(2, 8, 61, 36, 90, 10), the selected m0s are c(1, 2, 17, 18, 14, 2)
# Hamming dist: cluster sizes are c(92, 9, 42, 45, 19), the selected m0s are c(38, 2, 16, 8, 4)

## for clustering in groups
load('GDSC_clust_res_in_grps_jaccard.RData')
sample_size=100
iter = 200
count = 0
# post pred check
for (i in 1:sample_size) {
  for(j in 1:iter){
    real_log_lik = 0
    fake_log_lik = 0
    
    for(k in 1:6){
      cl = list(K = cl_samples[[k]]$KK[i], N = cl_samples[[k]]$NN[i, ], C = cl_samples[[k]]$CC[i, ])
      active_cid = which(col_block[k, ] > 0)
      cpds = colnames(col_block[, active_cid])
      A = GDSC[, cpds] 
      a = rep(mean(GDSC, na.rm = T), ncol(A))
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
print(count / sample_size / iter) ## 7 / 20000 for hamming dist; 2/20000 for jaccard dist, B = 20000
# avg log likelihood
real_log_lik = 0
for (i in 1:sample_size) {
  for(j in 1:iter){
    for(k in 1:5){
      cl = list(K = cl_samples[[k]]$KK[i], N = cl_samples[[k]]$NN[i, ], C = cl_samples[[k]]$CC[i, ])
      active_cid = which(col_block[k, ] > 0)
      cpds = colnames(col_block[, active_cid])
      A = GDSC[, cpds] 
      a = rep(mean(GDSC, na.rm = T), ncol(A))
      b = 1 - a
      post_phi = generate_posterior_parameters(cl, A, a, b)
      real_log_lik = real_log_lik + log_likelihood(cl, post_phi, mat = A, x0 = A)
    }
  }
}
print(real_log_lik / sample_size / iter) ## 7298.78 for jaccard spectral; 7261.77 for hamming spectral


# lof <- function(Z){
#   binary_str = apply(Z, 2, function(z){
#     return(paste(z, collapse = ""))
#   })
#   col_order = order(binary_str, decreasing = T)
#   return(Z[, col_order])
# }
# col_block = lof(col_block)
# k = 3
# block = GDSC[, names(col_block[k, col_block[k,]>0])]
# plot(block)

cl = clara(x = GDSC, k = 6, metric = 'jaccard', cluster.only = T)
cl = as.factor(cl)
summary(cl)
