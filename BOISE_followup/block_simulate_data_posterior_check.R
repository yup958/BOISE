setwd('~/RAwork/BOISE/BOISE_followup/')
### Define a function to find m0 that minimize prior and posterior difference
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
  min_diff = 1000
  best_m0 = 1
  
  for(m0 in seq(lower, upper, by = 3)){
    prior_m0_expect = m0 * harmonic_sum(start = m0, end = m0+n-1)
    cl_sample = dpmm_beta(x0, a, b, m0, burn_in = 500, sample_size=100, thinning = 5)
    diff = abs(mean(cl_sample$KK) - prior_m0_expect)
    if(diff < min_diff){
      min_diff = diff
      best_m0 = m0
    }
  }
  return(best_m0)
}

M_generate <- function(P){
  M = apply(P, 2, function(p){
    return(rbinom(length(p), 1, p))
  })
  M_sum = apply(M, 1, sum)
  while (length(which(M_sum == 0)) > 0) {
    M = apply(P, 2, function(p){
      return(rbinom(length(p), 1, p))
    })
    M_sum = apply(M, 1, sum)
  }
  return(M)
}

alpha = 0.06
beta = 0.94

set.seed(2)
# generating matrices with clustering in rows and columns
col_cl = BOISE::Initial_beta(x0 = matrix(0, 200, 200), m0 = 10)
P = matrix(0, nrow=100, ncol = 200)
prev_col = 0
for(k in 1:col_cl$K){
  ncols = col_cl$N[k]
  m0 = sample(3:15, size = 1, replace = T)
  row_cl = BOISE::Initial_beta(x0 = matrix(0, 100, 100), m0 = m0)
  probs = rbeta(row_cl$K, alpha, beta)
  prev_row = 0
  for(r in 1:row_cl$K){
    nrows = row_cl$N[r]
    P[(prev_row+1) : (prev_row+nrows), (prev_col+1):(prev_col+ncols)] = probs[r]
    prev_row = prev_row + nrows
  }
  prev_col = prev_col + ncols
}
plot(P, main = "Synthetic Probability Matrix with Blocks", xlab = 'Compounds', ylab = 'Targets')
block_x0 = M_generate(P)

### Posterior predictive check

### help functions
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

log_likelihood <- function(cl, Phi, mat){
  m = nrow(mat)
  neg_log_lik = 0
  for(i in 1:m){
    ck = cl$C[i]
    P = Phi[ck, ]
    neg_log_lik = neg_log_lik - sum(dbinom(mat[i,], 1, P, log = TRUE))
  }
  return(neg_log_lik)
}
a = rep(mean(block_x0),ncol(block_x0))
b = 1-a
sample_size = 1000

m0 = m0_Find(x0 = block_x0, lower = 3, upper = 15)
system.time({
  cl_sample = dpmm_beta(x0 = block_x0, alpha = a, beta = b, m0 = m0,
                        burn_in=500, sample_size = sample_size, thinning = 10)
})

count = 0
for (i in 1:sample_size) {
  cl = list(K = cl_sample$KK[i], N = cl_sample$NN[i, ], C = cl_sample$CC[i, ])
  post_phi = generate_posterior_parameters(cl, block_x0, a, b)
  Y = generate_fake_data(cl, post_phi)
  real_log_lik = log_likelihood(cl, post_phi, block_x0)
  fake_log_lik = log_likelihood(cl, post_phi, Y)
  if(fake_log_lik > real_log_lik){
    count = count + 1
  }
}
print(count / sample_size) ## 0.003, for col_cl m0 = 10 and row_cl m0 sample from 3 to 15.
